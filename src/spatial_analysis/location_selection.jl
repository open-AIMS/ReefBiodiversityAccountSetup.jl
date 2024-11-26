using DataFrames, GeoDataFrames, Distances

"""
    suggest_impact_sites(site_data::DataFrame; min_site_karea::Float64=1000.0, sorted=true)::DataFrame

Calculate measure of suitable sites for implementing activities given a set of spatial data

# Arguments
- `site_data` : Spatial dataset including geomorphic classes for each site.
- `min_site_karea` : Minimum acceptable k area for sites to implement on.
- `sorted` : true is output should be sorted from best to worst (don't use if using to plot map)
"""
function suggest_impact_sites(
    site_data::DataFrame;
    min_site_karea::Float64=1000.0,
    sorted=true
)::DataFrame
    n_locs = size(site_data, 1)
    locations = collect(1:n_locs)
    geomorphic_protection = zeros(n_locs)

    geom_vec = site_data[:, :class]
    geomorphic_protection[geom_vec .== "Sheltered Reef Slope"] .= 3
    geomorphic_protection[geom_vec .== "Back Reef Slope"] .= 2
    geomorphic_protection[geom_vec .== "Reef Slope"] .= 1
    geomorphic_protection[geom_vec .== "Deep Lagoon"] .= 1

    k_area = site_data[:, "k"] .* site_data[:, "area"]

    protection_rating =
        normalize(geomorphic_protection) + normalize(site_data[:, "depth_med"])
    protection_rating[k_area .< min_site_karea] .= 0.0

    rating = protection_rating ./ sum(protection_rating)
    if sorted
        s_order::Vector{Int64} = sortperm(rating; rev=true)

    else
        s_order = 1:n_locs
    end
    return DataFrame(
        hcat(locations[s_order], rating[s_order]), ["locations", "rating"])
end

"""
    suggest_control_sites(impact_site_id::Int64, site_data::DataFrame, category_constraints::Union{Vector{Symbol},Vector{String}}; weightings::Vector{Float64}=ones(size(site_data, 2) - (1 + length(category_constraints))), ID_COLUMN::Union{String,Symbol}=:reef_siteid, distance_func=chebyshev)::DataFrame

Output ordered list of sites which are similar to impact site. Similarity measured by
normalised Chebyshev distance (or `distance_func` where included).

# Arguments
- `impact_site_id` : Id of site to be used to implement activities at.
- `site_data` : Data to be used to judge similarity to impact site.
- `category_constraints` : List of column names in `site_data` which are categories that the
    impact and control sites must both sit within.
- `weightings` : Weightings for the similarity criteria, to weight importance in control site selection.
    Default is equal weighting (all 1.0).
- `ID_COLUMN` : Column in `site_data` which is the unique identifier for sites.
- `distance_func` : A distance function used to judge "similarity".
"""
function suggest_control_sites(
    impact_site_id::Int64,
    site_data::DataFrame,
    category_constraints::Union{Vector{Symbol},Vector{String}};
    weightings::Vector{Float64}=ones(
        size(site_data, 2) - (1 + length(category_constraints))
    ),
    ID_COLUMN::Union{String,Symbol}=:reef_siteid,
    distance_func=chebyshev)::DataFrame
    impact_site_data = site_data[impact_site_id, :]
    constraints::Vector{Bool} = fill(true, size(site_data, 1))
    for cc in category_constraints
        sites_subset = site_data[:, cc] .== impact_site_data[cc]
        constraints = constraints .& sites_subset
    end

    criteria_df = site_data[constraints, Not(ID_COLUMN, category_constraints...)]
    distances =
        distance_func.(
            Array(impact_site_data[Not(ID_COLUMN, category_constraints...)])',
            Matrix(criteria_df)
        )

    scores = normalize(distances) .* weightings
    distances = dropdims(sum(scores; dims=2); dims=2)
    s_order::Vector{Int64} = sortperm(distances; rev=false)

    return DataFrame(
        hcat(findall(constraints)[s_order], scores[s_order, :], distances[s_order]),
        vcat(["Location index"], names(criteria_df), ["Similarity"])
    )
end
