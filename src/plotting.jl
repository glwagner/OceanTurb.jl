module Plotting

export
    plot,
    removespine,
    removespines,
    cornerspines,
    bottomspine,

    defaultcolors

using OceanTurb, PyPlot, PyCall

# A few nice things for plotting
import PyPlot: plot

font_manager = pyimport("matplotlib.font_manager")
defaultcolors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

plot(f::AbstractField, args...; kwargs...) = plot(data(f), nodes(f), args...; kwargs...)
plot(op::Function, f::AbstractField, args...; kwargs...) = plot(op.(data(f)), nodes(f), args...; kwargs...)

"Remove `spine` from `ax`."
function removespine(side, ax=gca())
    ax.spines[side].set_visible(false)
    keywords = Dict(Symbol(side)=>false, Symbol(:label, side)=>false)
    ax.tick_params(keywords)
    nothing
end

removespines(sides...; ax=gca()) = for side in sides; removespine(side, ax); end
cornerspines(ax=gca()) = removespines("top", "right"; ax=ax)

function bottomspine(ax=gca())
    removespines("top", "right", "left"; ax=ax)
    ax.tick_params(left=false, labelleft=false)
end

end # module
