import GeoDataFrames as GDF
using DataFrames, YAXArrays
using CSV, MAT
using NCDatasets
import ArchGDAL as AG
using Glob
using GeoFormatTypes: EPSG
using Rasters
using Rasters: Intervals
using Statistics
using TOML

function load_config(; config_path::String="config.toml")::Dict
    return TOML.parsefile(config_path)
end

"""
    load_spatial_base(config_file::Dict)

Loads key spatial data layers from an Allen Atlas dataset (benthic, geomorphic and extent (experimental)).

# Returns
- `benthic` : Dataframe of polygons and their benthic classes
- `geomorphic` : Dataframe of polygons and their geomorphic classes
- `reef_extent` : Dataframe of polygons showing reef extent (experiemental product in the Allen Atlas)
"""
function load_spatial_base(config_file::Dict)
    # Base allen coral atlas data
    allen_data_base = config_file["spatial_data"]["allen_spatial"]

    # Load benthic data
    benthic = GDF.read(string(allen_data_base, "\\Benthic-Map\\benthic.gpkg"))

    # Load geomorphic data
    geomorphic = GDF.read(
        string(allen_data_base, "\\Geomorphic-Map\\geomorphic.gpkg")
    )
    # Load reef extent data
    reef_extent = GDF.read(
        string(allen_data_base, "\\Reef-Extent\\reefextent.gpkg")
    )

    return benthic, geomorphic, reef_extent
end

"""
    get_geo_within_box(geo_data::DataFrame, box_upper::Tuple, box_lower::Tuple)::DataFrame

Get `geo_data` within specified bounding box.

# Arguments
- `geo_data` : GeoDataFrame
- `box_upper` : Upper coordinates for bounding box
- `box_lower` : Lowert coordinates for bounding box
"""
function get_geo_within_box(
    geo_data::DataFrame,
    box_upper::Tuple,
    box_lower::Tuple
)::DataFrame
    centroids = AG.centroid(geo_data[:, :geom])
    longs = AG.getx.(centroids, 0)
    lats = AG.gety.(centroids, 0)
    in_box =
        (lats .<= box_upper[1]) .& (lats .>= box_lower[1]) .& (longs .>= box_upper[2]) .&
        (longs .<= box_lower[2])
    return geo_data[in_box, :]
end

"""
    filter_site_area(df::DataFrame)::DataFrame

Filter benthic data to only include polygons suitable for coral growth.

# Arguments
- `df` : spatial data include benthic class column
"""
function filter_site_area(df::DataFrame)::DataFrame
    return df[(df.class .== "Coral/Algae") .|| (df.class .== "Rock"), :]
end

"""
    get_geomorphology_coral_area(df_1::DataFrame, df_2::DataFrame)::DataFrame

Equivalent to a spatial join for df_1 and df_2, where the operation is intersection.

# Arguments
- `df_1` : Inner dataframe
- `df_2` : Outer dataframe
"""
function get_geomorphology_coral_area(
    df_1::DataFrame,
    df_2::DataFrame
)::DataFrame
    geo_join = crossjoin(df_1, df_2; makeunique=true)
    subset!(geo_join, [:geom, :geom_1] => (a, b) -> intersects.(a, b))
    return select!(geo_join, Not(:geom_1))
end

"""
    multipoly_geom_within(outer_gdf::DataFrame, inner_gdf::DataFrame, classes::Union{String,Symbol})

Combine multipolygons in `outer_gdf` of a particular class which sit inside the polygons of `inner_gdf`

# Arguments
- `outer_gdf` : Outer dataframe
- `inner_gdf` : Inner dataframe
- `classes` : Vector giving the unique class of each polygon in `inner_gdf`
"""
function multipoly_geom_within(
    outer_gdf::DataFrame,
    inner_gdf::DataFrame,
    classes::Union{String,Symbol}
)
    new_df_store = DataFrame(; geom=AG.createmultipolygon(), class=String[])
    # add to this to get multipolys for each geomorphic class
    for poly in outer_gdf[:, :geom]
        poly_temp = DataFrames.groupby(
            inner_gdf[findall(AG.within.(inner_gdf[:, :geom], [poly])), [:geom, classes]],
            classes
        )
        if !isempty(poly_temp)
            for df in poly_temp
                multipoly_temp = AG.createmultipolygon()
                for pp in df[:, :geom]
                    AG.addgeom!(multipoly_temp, pp)
                end
                push!(new_df_store, [multipoly_temp, df[:, classes][1]])
            end
        end
    end
    return new_df_store
end

"""
    multipoly_geom_intersection(outer_gdf::DataFrame, inner_gdf::DataFrame, classes::Union{String,Symbol})

Find the intersections of `outer_gdf` with each poly in `inner_gdf` and combine to create new set of geometries
with unique classes in `inner_gdf`.

# Arguments
- `outer_gdf` : Outer dataframe
- `inner_gdf` : Inner dataframe
- `classes` : Vector giving the unique class of each polygon in `inner_gdf`
"""
function multipoly_geom_intersection(
    outer_gdf::DataFrame,
    inner_gdf::DataFrame,
    classes::Union{String,Symbol}
)
    new_df_store = DataFrame(; geom=AG.createmultipolygon(), class=String[])
    # add to this to get multipolys for each geomorphic class
    for poly in outer_gdf[:, :geom]
        intersect_idx = findall(AG.intersects.(
            inner_gdf[:, :geom], [poly]
        ))

        if !isempty(intersect_idx)
            poly_temp = AG.intersection.(inner_gdf[:, :geom], [poly])[intersect_idx]
            classes_temp = inner_gdf[intersect_idx, classes]
            for class in unique(classes_temp)
                multipoly_temp = AG.createmultipolygon()
                for pp in poly_temp[classes_temp .== class]
                    if (typeof(pp) == AG.IGeometry{AG.wkbGeometryCollection}) ||
                        (typeof(pp) == AG.IGeometry{AG.wkbMultiPolygon})
                        for p_sub in 0:(AG.ngeom(pp) - 1)
                            if typeof(AG.getgeom(pp, p_sub)) == AG.IGeometry{AG.wkbPolygon}
                                AG.addgeom!(multipoly_temp, AG.getgeom(pp, p_sub))
                            end
                        end
                    elseif (typeof(pp) == AG.IGeometry{AG.wkbPolygon})
                        AG.addgeom!(multipoly_temp, pp)
                    end
                end
                if !AG.isempty(multipoly_temp)
                    push!(new_df_store, [multipoly_temp, class])
                end
            end
        end
    end
    return new_df_store
end

"""
    get_multipoly_area(outer_gdf::DataFrame, inner_gdf::DataFrame)

Find the area of `inner_gdf` intersecting with `outer_gdf`.

# Arguments
- `outer_gdf` : Outer dataframe
- `inner_gdf` : Inner dataframe
"""
function get_multipoly_area(
    outer_gdf::DataFrame,
    inner_gdf::DataFrame
)
    multipoly_area_store = zeros(Float64, length(outer_gdf[:, :geom]))
    proj_str = ProjString("+proj=utm +zone=55 +south +datum=WGS84 +units=m +no_defs")
    for (ind_p, poly) in enumerate(outer_gdf[:, :geom])
        temp_poly = AG.intersection.(inner_gdf[:, :geom], [poly])[findall(
            AG.intersects.(inner_gdf[:, :geom], [poly])
        )]
        if !isempty(temp_poly)
            proj_m = AG.reproject(temp_poly, EPSG(4326), proj_str; order=:trad)
            multipoly_area_store[ind_p] = sum(AG.geomarea.(proj_m))
        end
    end
    return multipoly_area_store
end

"""
    set_reef_k(geomorphic_gdf::DataFrame, benthic_gdf::DataFrame; storage_string::String="spatial_data_temp.gpkg", threshold_k::Float64=0.05)

Set the proportional k area and total site area for the set of polygons in `reef_gdf`, given a set of k areas
    (calculated using `get_multipoly_area`).

# Arguments
- `geomorphic_gdf` : geomorphic classes and polygons.
- `benthic_gdf` : benthic classes and polygons.
- `storage_string` : filename for temporary file saving to allow coordinate transforms and other manipulations.
- `threshold_k` : threshold for k below which polygons are filtered out
"""
function set_reef_k(
    geomorphic_gdf::DataFrame,
    benthic_gdf::DataFrame;
    storage_string::String="spatial_data_temp.gpkg",
    threshold_k::Float64=0.05
)::DataFrame
    benthic_filtered = benthic_gdf[
        benthic_gdf[:, "class"] .== "Coral/Algae", :
    ]

    # Get area available for Coral/Algal type on polygons
    areas_k = get_multipoly_area(geomorphic_gdf, benthic_filtered)

    # Save gdf as temporary file to allow projection
    GDF.write(storage_string, geomorphic_gdf; driver="GPKG", geom_columns=(:geom,))
    gdf_temp = GDF.read(storage_string)

    proj_str = ProjString("+proj=utm +zone=55 +south +datum=WGS84 +units=m +no_defs")
    proj_m = AG.reproject(gdf_temp[:, :geom], EPSG(4326), proj_str; order=:trad)

    areas = AG.geomarea.(proj_m)
    k = areas_k ./ areas
    geomorphic_gdf[!, :area] = areas
    geomorphic_gdf[!, :k] = k
    return geomorphic_gdf[geomorphic_gdf.k.>threshold_k, :]
end

"""
    get_depths(reef_gdf::DataFrame, config_file::Dict; temporary_gpkg_name::String="spatial_data_temp.gpkg",)::Tuple{DataFrame,Array}
    get_depths(reef_gdf::DataFrame, gbrmpa_region_path::String, gbrmpa_bathy_path::String; temporary_gpkg_name::String="spatial_data_temp.gpkg",)::Tuple{DataFrame,Array}

Get depths using data sourced from GBRMPA.

# Arguments
- `config_file` : Dict of key filenames loaded from a "config.toml file"
- `bathy_data_dir` : Directory for GBRMPA bathymetry data
- `temporary_gpkg_name` : Filename for temporary geopackage file saved in any previous steps.
- `region_path` : Location of file containing regional GBRMPA spatial data.
"""
function get_depths(
    reef_gdf::DataFrame,
    config_file::Dict;
    temporary_gpkg_name::String="spatial_data_temp.gpkg"
)::Tuple{DataFrame,Array}
    return get_depths(
        reef_gdf,
        config_file["spatial_data"]["gbrmpa_region"],
        config_file["bathy_data"]["gbrmpa_bathy"];
        temporary_gpkg_name=temporary_gpkg_name
    )
end
function get_depths(
    reef_gdf::DataFrame,
    gbrmpa_region_path::String,
    gbrmpa_bathy_path::String;
    temporary_gpkg_name::String="spatial_data_temp.gpkg"
)::Tuple{DataFrame,Array}
    function summary_func(x)
        return if (length(unique(x)) > 1)
            [-maximum(x), -mean(x), -median(x), -minimum(x), std(x)]
        else
            first(x)
        end
    end
    # Save gdf as temporary file to allow projection
    GDF.write(temporary_gpkg_name, reef_gdf; driver="GPKG", geom_columns=(:geom,))
    gdf = GDF.read(temporary_gpkg_name)
    depths = zeros(Float64, size(gdf, 1), 5)
    errored_empty = zeros(Int64, size(gdf, 1))

    REGIONS = String[
        "Townsville-Whitsunday",
        "Cairns-Cooktown",
        "Mackay-Capricorn",
        "FarNorthern"
    ]

    for reg in REGIONS
        src_bathy_path = first(glob("*.tif", joinpath(gbrmpa_bathy_path, reg)))
        src_bathy = Raster(src_bathy_path; mappedcrs=EPSG(4326), lazy=true)

        proj_str = ProjString(
            AG.toPROJ4(AG.importWKT(crs(src_bathy).val; order=:compliant))
        )

        # Ensure polygon types match
        region_features = GDF.read(gbrmpa_region_path)
        # Force CRS to match raster data
        region_features.geometry = AG.reproject(
            region_features.SHAPE, EPSG(4326), proj_str; order=:trad
        )
        region_features[!, :geometry] = Vector{AG.IGeometry}(
            AG.forceto.(region_features.geometry, AG.wkbMultiPolygon)
        )

        reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)
        tgt_region = region_features.geometry[reg_idx]

        gdf = GDF.read(temporary_gpkg_name)
        gdf[!, :geom] = Vector{AG.IGeometry}(AG.forceto.(gdf.geom, AG.wkbMultiPolygon))
        target_geoms = gdf.geom

        target_geoms = AG.reproject(target_geoms, EPSG(4326), proj_str; order=:trad)

        feature_match_ids = unique(findall(GDF.intersects.(tgt_region, target_geoms)))

        if !isempty(feature_match_ids)
            for id in feature_match_ids
                try
                    depths[id, :] .= zonal(summary_func, src_bathy; of=target_geoms[id])
                catch err
                    if (err isa MethodError) || (err isa ArgumentError)
                        # Raises MethodError when `target_geoms` is empty
                        # Raises ArgumentError where `zonal()` produces `Missing` (no data)
                        msg = "MethodError or ArgumentError on $(id)\n"
                        msg = msg * "Possibly a reef feature with no overlapping polygon"
                        errored_empty[id] = 1
                        continue
                    end

                    @info "Error on $(id)"
                    rethrow(err)
                end
            end
        end
    end
    depths[depths .< 0.0] .= 7.0
    reef_gdf[!, :depth_med] = depths[:, 3]

    return reef_gdf, depths
end

"""
    median_features_allen(reef_gdf::DataFrame, config_file::Dict; temporary_gpkg_name::String="spatial_data_temp.gpkg", data_name::Symbol=:depth_med, is_depth=false)::Tuple{DataFrame,Array}
    median_features_allen(reef_gdf::DataFrame, allen_bathy_filepath::String; temporary_gpkg_name::String="spatial_data_temp.gpkg", data_name::Symbol=:depth_med, is_depth=false)::Tuple{DataFrame,Array}

Get median values of an Allen Atlas Raster over a set of geometries

# Arguments
- `config_file` : Dict of key filenames loaded from a "config.toml file"
- `reef_gdf` : Current reef geodataframe to augment
- `allen_dir` : Directory for Allen Atlas Raster data (net cdf or tif)
- `temporary_gpkg_name` : Filename for temporary geopackage file saved in any previous steps
"""
function median_features_allen(
    reef_gdf::DataFrame,
    config_file::Dict;
    temporary_gpkg_name::String="spatial_data_temp.gpkg",
    data_name::Symbol=:depth_med,
    is_depth=false
)::Tuple{DataFrame,Array}
    return median_features_allen(
        reef_gdf,
        config_file["bathy_data"]["allen_bathy"];
        temporary_gpkg_name=temporary_gpkg_name,
        data_name=data_name,
        is_depth=is_depth
    )
end
function median_features_allen(
    reef_gdf::DataFrame,
    allen_bathy_filepath::String;
    temporary_gpkg_name::String="spatial_data_temp.gpkg",
    data_name::Symbol=:depth_med,
    is_depth=false
)::Tuple{DataFrame,Array}
    # Save temporary file to allow projections to retrieve depths
    GDF.write(temporary_gpkg_name, reef_gdf; driver="GPKG", geom_columns=(:geom,))

    src_allen = Raster(allen_bathy_filepath; mappedcrs=EPSG(4326), lazy=true)
    proj_str = ProjString(AG.toPROJ4(AG.importWKT(crs(src_allen).val; order=:compliant)))

    function summary_func(x)
        return if (length(unique(x)) > 1)
            [maximum(x), mean(x), median(x), minimum(x), std(x)]
        else
            first(x)
        end
    end

    gdf = GDF.read(temporary_gpkg_name)
    gdf[!, :geom] = Vector{AG.IGeometry}(AG.forceto.(gdf.geom, AG.wkbMultiPolygon))
    target_geoms = gdf.geom
    target_geoms = AG.reproject(target_geoms, EPSG(4326), proj_str; order=:trad)

    data_store = zeros(Float64, size(gdf, 1), 5)
    errored_empty = zeros(Int64, size(gdf, 1))

    for (id, id_geom) in enumerate(target_geoms)
        try
            data_store[id, :] .= zonal(summary_func, src_allen; of=id_geom, shape=:polygon)
        catch err
            if (err isa MethodError) || (err isa ArgumentError)
                # Raises MethodError when `target_geoms` is empty
                # Raises ArgumentError where `zonal()` produces `Missing` (no data)
                msg = "MethodError or ArgumentError on $(id)\n"
                msg = msg * "Possibly a reef feature with no overlapping polygon"
                errored_empty[id] = 1
                continue
            end

            @info "Error on $(id)"
            rethrow(err)
        end
    end
    if is_depth
        data_store = data_store ./ 100 # convert to m from cm
        data_store[data_store .<= 0.0] .= 7.0
    end
    reef_gdf[!, data_name] = data_store[:, 3]

    return reef_gdf, data_store
end

"""
    noaa_dhw_means(reef_gdf::DataFrame, config_file::Dict; temporary_gpkg_name::String="geomorph_temp.gpkg")::Tuple{DataFrame,YAXArray}
    noaa_dhw_means(reef_gdf::DataFrame, dhw_filepath::String; temporary_gpkg_name::String="geomorph_temp.gpkg")::Tuple{DataFrame,YAXArray}

Get zonal mean values of the NOAA DHW product over a set of geometries in a geopackage.

# Arguments
- `config_file` : Dict of key filenames loaded from a "config.toml file"
- `reef_gdf` : Reef geodataframe to augment
- `dhw_fn` : Directory for NOAA DHW product data
- `temporary_gpkg_name` : Filename for temporary geopackage file saved in any previous steps.
"""
function noaa_dhw_means(
    reef_gdf::DataFrame,
    config_file::Dict;
    temporary_gpkg_name::String="geomorph_temp.gpkg"
)::Tuple{DataFrame,YAXArray}
    return noaa_dhw_means(
        reef_gdf,
        config_file["dhw_data"]["noaa_dhw"];
        temporary_gpkg_name=temporary_gpkg_name
    )
end
function noaa_dhw_means(
    reef_gdf::DataFrame,
    dhw_filepath::String;
    temporary_gpkg_name::String="geomorph_temp.gpkg"
)::Tuple{DataFrame,YAXArray}
    GDF.write(temporary_gpkg_name, reef_gdf; driver="GPKG", geom_columns=(:geom,))

    noaa_raster = Raster(
        dhw_filepath;
        mappedcrs=EPSG(4326),
        crs=convert(WellKnownText, EPSG(4326)),
        name="ann_max_dhw",
        lazy=true
    )
    noaa_raster = set(noaa_raster, X => Intervals)
    noaa_raster = set(noaa_raster, Y => Intervals)

    proj_str = ProjString(AG.toPROJ4(AG.importWKT(crs(noaa_raster).val; order=:compliant)))
    gdf = GDF.read(temporary_gpkg_name)
    gdf[!, :geom] = Vector{AG.IGeometry}(AG.forceto.(gdf.geom, AG.wkbMultiPolygon))
    target_geoms = gdf.geom
    target_geoms = AG.reproject(target_geoms, EPSG(4326), proj_str; order=:trad)

    dhw_means = zeros(Float64, length(noaa_raster.dims[3]), size(gdf, 1))
    errored_empty = zeros(Int64, length(noaa_raster.dims[3]))

    for (yr_id, yr) in enumerate(collect(noaa_raster.dims[3]))
        try
            temp_raster = resample(noaa_raster[Ti=At(yr)]; crs=crs(target_geoms[1]))
            dhw_means[yr_id, :] .= zonal(
                mean, temp_raster; of=target_geoms, shape=:polygon, boundary=:touches
            )
        catch err
            if (err isa MethodError) || (err isa ArgumentError)
                # Raises MethodError when `target_geoms` is empty
                # Raises ArgumentError where `zonal()` produces `Missing` (no data)
                msg = "MethodError or ArgumentError on $(yr_id)\n"
                msg = msg * "Possibly a reef feature with no overlapping polygon"
                errored_empty[yr_id] = 1
                continue
            end

            @info "Error on $(yr_id)"
            rethrow(err)
        end
    end

    ax_list = (
        Dim{:timesteps}(noaa_raster.dims[3].val.data), Dim{:sites}(1:length(target_geoms))
    )
    reef_gdf[!, :dhw_hist_mean] = vec(dropdims(mean(dhw_means; dims=1); dims=1))
    reef_gdf[!, :dhw_hist_sd] = vec(dropdims(std(dhw_means; dims=1); dims=1))

    return reef_gdf, YAXArray(ax_list, dhw_means)
end

"""
    normalize(x)

Min-max normalisation of a vector
"""
function normalize(x)
    return (x .- minimum(x)) ./ (maximum(x) - minimum(x))
end
