using GeoDataFrames
import ArchGDAL as AG
using DataFrames, YAXArrays
using Statistics
using GLMakie, GeoMakie, GraphMakie
using DimensionalData

"""
    spatial_map(geo_df::DataFrame, color_vec::Vector{String}; fig_opts::Dict{Symbol,Any}=set_figure_defaults(Dict()), axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()), opts::Dict{Symbol,Any}=Dict())::Figure
    spatial_map(f::Figure, geo_df::DataFrame, color_vec::Vector{String}; axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()), opts::Dict{Symbol,Any}=Dict())::Figure
    spatial_map(geo_df::DataFrame, color_vec::Union{Vector{Float32},Vector{Float64}}; axis_opts=Dict{Symbol,Any}=set_axis_defaults(Dict()), opts::Dict{Symbol,Any}=Dict())::Figure
    spatial_map(f::Figure, geo_df::DataFrame, color_vec::Union{Vector{Float32},Vector{Float64}}; axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()), opts::Dict{Symbol,Any}=Dict())::Figure

Plots a spatial heat map using a color vector of floats (for a continuous mapping) or a color vector of strings (for categorical mapping)

# Arguments
- `geo_df` : GeoDataFrame containing geometries to be plotted
- `color_vec` : Vector of values of length of `geo_df`
- `opts` : Aviz options
    - `legend_name` : name for legend showing categories and their colors
    - `colorbar_label` : name for colorbar showing color range
    - `color_map` : discretized colormap for categories or continuous for heatmaps
    - `color_range` : Numerical limits for continuous colorbar
- `axis_opts` : Additional options to pass to adjust Axis attributes
  See: https://docs.makie.org/v0.19/api/index.html#Axis
"""
function spatial_map(
    geo_df::DataFrame,
    color_vec::Vector{String};
    fig_opts::Dict{Symbol,Any}=set_figure_defaults(Dict()),
    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    f = Figure(; fig_opts...)
    return spatial_map(f, geo_df, color_vec; axis_opts=axis_opts, opts=opts)
end
function spatial_map(
    f::Figure,
    geo_df::DataFrame,
    color_vec::Vector{String};
    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    legend_name = get(opts, :legend_name, "Class")
    color_map = get(opts, :color_map, :batlow25)

    levels = unique(color_vec)

    colors = Makie.categorical_colors(color_map, length(levels))
    legend_elems::Vector{Any} = ones(length(levels))
    spatial = GeoAxis(
        f[1, 1];
        axis_opts...
    )

    for (int_lev, lev) in enumerate(levels)
        temp_geo = geo_df[findall(color_vec .== lev), :]
        geo_data = GeoMakie.to_multipoly(temp_geo[:, :geom])

        poly!(
            spatial,
            geo_data;
            color=colors[int_lev],
            label=lev,
            transparency=true
        )
    end
    legend_elems = [PolyElement(; color=colors[int_lev]) for int_lev in 1:length(levels)]
    Legend(f[1, 2], legend_elems, levels, legend_name)

    return f
end
function spatial_map(
    geo_df::DataFrame,
    color_vec::Union{Vector{Float32},Vector{Float64}};
    fig_opts::Dict{Symbol,Any}=set_figure_defaults(Dict()),
    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    f = Figure(; fig_opts...)
    return spatial_map(f, geo_df, color_vec; axis_opts=axis_opts, opts=opts)
end
function spatial_map(
    f::Figure,
    geo_df::DataFrame,
    color_vec::Union{Vector{Float32},Vector{Float64}};
    axis_opts::Dict{Symbol,Any}=set_axis_defaults(Dict()),
    opts::Dict{Symbol,<:Any}=Dict{Symbol,Any}()
)::Figure
    colorbar_label = get(opts, :colorbar_label, "Metric")
    color_map = get(opts, :color_map, :batlow25)
    color_range = get(opts, :color_range, (0.0, maximum(color_vec)))

    geo_data = GeoMakie.to_multipoly(geo_df[:, :geom])
    spatial = GeoAxis(
        f[1, 1];
        axis_opts...
    )
    poly!(
        spatial,
        geo_data;
        color=color_vec,
        colormap=color_map,
        colorrange=color_range
    )
    Colorbar(
        f[1, 2];
        colorrange=color_range,
        colormap=color_map,
        label=colorbar_label,
        height=Relative(0.70)
    )
    return f
end

"""
    set_figure_defaults(fig_opts::Dict{Any,Any})::Dict{Symbol,Any}

Set default figure settings for spatial figures
"""
function set_figure_defaults(
    fig_opts::Dict{Any,Any}
)::Dict{Symbol,Any}
    fig_opts[:size] = get(fig_opts, :size, (600, 900))
    fig_opts[:xticklabelsize] = get(fig_opts, :xticklabelsize, 14)
    fig_opts[:yticklabelsize] = get(fig_opts, :yticklabelsize, 14)

    return fig_opts
end

"""
    set_axis_defaults(axis_opts::Dict{Any,Any})::Dict{Symbol,Any}

Set default axis settings for spatial figures
"""
function set_axis_defaults(
    axis_opts::Dict{Any,Any}
)::Dict{Symbol,Any}
    axis_opts[:title] = get(axis_opts, :title, "Study Area")
    axis_opts[:xlabel] = get(axis_opts, :xlabel, "Longitude")
    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Latitude")
    axis_opts[:xgridwidth] = get(axis_opts, :xgridwidth, 0.5)
    axis_opts[:ygridwidth] = get(axis_opts, :ygridwidth, 0.5)
    axis_opts[:dest] = get(axis_opts, :dest, "+proj=latlong +datum=WGS84")

    return axis_opts
end
