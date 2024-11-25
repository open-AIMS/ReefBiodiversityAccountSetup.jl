module ReefBiodiversityAccountSetup

export suggest_impact_sites, suggest_control_sites
export load_config, load_spatial_base
export get_geo_within_box,
    get_multipoly_geom, get_multipoly_geom_intersection, get_multipoly_area
export filter_site_area, get_geomorphology_coral_area, set_reef_k, get_depths,
    get_median_features_allen, noaa_dhw_means

# Spatial analysis
include("./spatial_analysis/spatial_analysis.jl")

export spatial_map, temporal_spread

# Account plotting
include("./plotting/plotting.jl")

export spatial_analysis, plotting

end
