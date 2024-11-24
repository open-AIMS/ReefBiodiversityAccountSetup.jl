using DataFrames

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

    new_df = DataFrame(
        hcat(locations, protection_rating ./ sum(protection_rating)),
        ["locations", "rating"]
    )
    if sorted
        return sort!(new_df, ["rating"]; rev=true)
    else
        return new_df
    end
end

"""
    suggest_control_sites(impact_site_id::Int64, site_data::DataFrame)::Matrix

Output ordered list of sites which are similar to impact site. Similarity measured by
normalised Chebyshev distance.

# Arguments
- `impact_site_id` : Id of site to be used to implement activities at.
- `site_data` : Data to be used to judge similarity to impact site. Must include geomorphic class.
"""
function suggest_control_sites(
    impact_site_id::Int64,
    site_data::DataFrame
)::Matrix
    impact_site_data = site_data[impact_site_id, Not(:geom)]
    geomorphic_class = impact_site_data.class
    sites_subset = findall(site_data.class .== geomorphic_class)

    distances =
        chebyshev.(
            Array(impact_site_data[Not(:class)])',
            Matrix(site_data[sites_subset, Not(:class, :geom)])
        )
    distances = dropdims(sum(normalize(distances); dims=2); dims=2)

    s_order::Vector{Int64} = sortperm(distances; rev=false)

    return Union{Float64,Int64}[Int64.(sites_subset[s_order]) distances[s_order]]
end
