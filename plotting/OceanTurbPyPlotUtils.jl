module OceanTurbPyPlotUtils

export
    plot,
    removespine,
    removespines,
    cornerspines,
    bottomspine,

    defaultcolors,
    fontmanager,
    usecmbright

using OceanTurb, PyPlot, PyCall

# A few nice things for plotting
import PyPlot: plot

fontmanager = pyimport("matplotlib.font_manager")
defaultcolors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

plot(f::AbstractField, args...; kwargs...) = plot(data(f), nodes(f), args...; kwargs...)
plot(op::Function, f::AbstractField, args...; kwargs...) = plot(op.(data(f)), nodes(f), args...; kwargs...)

latexpreamble = """
\\usepackage{cmbright}
\\renewcommand{\\b}[1]    {\\boldsymbol{#1}}
\\renewcommand{\\r}[1]    {\\mathrm{#1}}
\\renewcommand{\\d}       {\\partial}
"""

function usecmbright()
  rc("text.latex", preamble=latexpreamble)
  rc("font", family="sans-serif")
  nothing
end

"Remove `spine` from `ax`."
removespine(side, ax=gca()) = ax.spines[side].set_visible(false)
removespines(sides...; ax=gca()) = for side in sides; removespine(side, ax); end

cornerspines(ax=gca()) = removespines("top", "right"; ax=ax)
bottomspine(ax=gca()) = removespines("top", "right", "left"; ax=ax)

end # module
