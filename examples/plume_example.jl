using OceanTurb, Printf

using OceanTurb.ModularKPP: DiagnosticPlumeModel

@use_pyplot_utils

usecmbright()

modelsetup = (N=16, L=128, stepper=:BackwardEuler)

Fb = 1e-7
Fu = 0.0
N² = 1e-5
Δt = 2minute
t_plot = (1hour, 2hour, 3hour)

name = "Free convection
    \\small{with \$ \\overline{w b} |_{z=0} = 10^{-7} \\, \\mathrm{m^2 \\, s^{-3}}\$}"

# The standard setup except with a plume model rather than a counter-gradient flux model.
model = ModularKPP.Model(; N=1024, L=64, stepper=:BackwardEuler, 
                         nonlocalflux = DiagnosticPlumeModel())

# Initial condition and fluxes
dTdz = N² / (model.constants.α * model.constants.g)
T₀(z) = 20 + dTdz * z

model.solution.T = T₀

Fθ = Fb / (model.constants.α * model.constants.g)
model.bcs.U.top = FluxBoundaryCondition(Fu)

model.bcs.T.top = FluxBoundaryCondition(Fθ)
model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

fig, axs = subplots(ncols=3)

sca(axs[1])
removespines("top", "right")
xlabel("Temperature \$ \\, {}^\\circ \\mathrm{C} \$")
ylabel(L"z \, \mathrm{(m)}")

sca(axs[3])
removespines("top", "right", "left")
axs[3].tick_params(left=false, labelleft=false)
xlabel("Vertical momentum")

sca(axs[2])
removespines("top", "right", "left")
axs[2].tick_params(left=false, labelleft=false)
xlabel("Plume temperature excess")

function iterate_and_plot(Δt, n=1)
    iterate!(model, Δt=Δt, Nt=1)

    sca(axs[1])
    plot(model.solution.T, "-", alpha=0.8, markersize=1.5)
    plot(model.state.plume.T, "--", alpha=0.8, markersize=1.5)

    ΔT = model.state.plume.T - model.solution.T
    sca(axs[2])
    plot(ΔT, "-", alpha=0.8, markersize=1.5)

    sca(axs[3])
    plot(model.state.plume.W², "-")

    return nothing
end

