using DataFrames, YAXArrays
using Statistics
using GLMakie, GeoMakie, GraphMakie
using DimensionalData

"""
    temporal_spread(f::Figure, data::YAXArray; axis_opts::Dict{Symbol,<:Any}=set_figure_defaults(Dict{Symbol,Any}()), opts::Dict{Symbol,<:Any}=set_figure_defaults(Dict{Symbol,Any}()))

Plot median of a data set over time with colored bands representing the 1st, 2nd and 3rd quantiles.

# Arguments
- `data` : Must include a :timesteps dimension + 2 other dims
- `opts` : Plotting options
    -`plot_color` : Color of lines and bands
- `axis_opts` : Additional options to pass to adjust Axis attributes
    -`ylabel` : Label for the metric being plotted
  See: https://docs.makie.org/v0.19/api/index.html#Axis
- `series_opts` : Additional options to pass to adjust Series attributes
  See: https://docs.makie.org/v0.19/api/index.html#series!
"""
function temporal_spread(
    f::Figure, data::YAXArray;
    axis_opts::Dict{Symbol,<:Any}=set_figure_defaults(Dict{Symbol,Any}()),
    opts::Dict{Symbol,<:Any}=set_figure_defaults(Dict{Symbol,Any}())
)
    n_timesteps = length(data.timesteps)
    data_quants = zeros(n_timesteps, 7)

    for tt ∈ 1:n_timesteps
        data_quants[tt, :] = quantile(
            data[tt, :, :], [0.15, 0.25, 0.35, 0.5, 0.65, 0.75, 0.85]
        )
    end

    axis_opts[:ylabel] = get(axis_opts, :ylabel, "Metric")
    plot_color = get(opts, :plot_color, :dodgerblue)
    years = collect(2025:2099)

    ax = Axis(
        f[1, 1];
        xlabel="Year",
        xlabelsize=20,
        ylabelsize=20,
        axis_opts...
    )

    band!(
        ax,
        years,
        vec(data_quants[:, 1]),
        vec(data_quants[:, 7]);
        color=plot_color,
        alpha=0.1
    )

    band!(
        ax,
        years,
        vec(data_quants[:, 2]),
        vec(data_quants[:, 6]);
        color=plot_color,
        alpha=0.3
    )

    band!(
        ax,
        years,
        vec(data_quants[:, 3]),
        vec(data_quants[:, 5]);
        color=plot_color,
        alpha=0.4
    )

    lines!(ax,
        years,
        vec(data_quants[:, 4]);
        color=plot_color,
        linewidth=5
    )
    return f
end
