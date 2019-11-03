using OceanTurb, Printf

using OceanTurb.ModularKPP: DiagnosticPlumeModel

@use_pyplot_utils

usecmbright()

modelsetup = (N=256, L=256, stepper=:BackwardEuler)

Fb = 1e-7
Fu = 0.0
N² = 1e-5
Δt = 2minute
t_plot = (0hour, 12hour, 48hour)

name = "Free convection
    \\small{with \$ \\overline{w b} |_{z=0} = 10^{-7} \\, \\mathrm{m^2 \\, s^{-3}}\$}"

# The standard setup except with a plume model rather than a counter-gradient flux model.
model = ModularKPP.Model(; modelsetup...,
                           nonlocalflux = DiagnosticPlumeModel())

# Initial condition and fluxes
dTdz = N² / (model.constants.α * model.constants.g)
T₀(z) = 20 + dTdz*z

model.solution.T = T₀

Fθ = Fb / (model.constants.α * model.constants.g)
model.bcs.U.top = FluxBoundaryCondition(Fu)

model.bcs.T.top = FluxBoundaryCondition(Fθ)
model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

iterate!(model, Δt)

#ModularKPP.update_state!(model)

ΔT = model.state.plumes.T - model.solution.T

fig, axs = subplots(ncols=3)

sca(axs[1])
plot(ΔT)

sca(axs[2])
plot(model.solution.T)
plot(model.state.plumes.T)

sca(axs[3])
plot(model.state.plumes.W²)
