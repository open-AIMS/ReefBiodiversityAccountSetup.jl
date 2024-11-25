module spatial_analysis

# Functions to import and process spatial layers from the Allen Atlas, GBRMPA and NOAA
include("process_layers.jl")
# Functions for selecting impact and control sites based on criteria
include("location_selection.jl")

end