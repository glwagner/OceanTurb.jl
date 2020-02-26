using OceanTurb, Printf 

using OceanTurb.ModularKPP: DiagnosticPlumeModel

@use_pyplot_utils

usecmbright()

function makeplot!(axs, model)
    ΔT = model.state.plume.T - model.solution.T

    ϵ = CellField(model.grid)
    for i in eachindex(ϵ)
        @inbounds ϵ[i] = ModularKPP.entrainment(model.grid.zc[i], model)
    end

    for ax in axs
        sca(ax); cla()
    end

    sca(axs[1])
    plot(ΔT)

    sca(axs[2])
    plot(model.solution.T)
    plot(model.state.plume.T, "--")

    sca(axs[3])
    plot(model.state.plume.W²)

    sca(axs[4])
    plot(ϵ)

    return nothing
end

modelsetup = (N=256, L=256, stepper=:BackwardEuler)

Fb = 1e-7
Fu = 0.0
N² = 1e-5
Δt = 1minute
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
OceanTurb.set!(model.state.plume.T, T₀)

Fθ = Fb / (model.constants.α * model.constants.g)
model.bcs.U.top = FluxBoundaryCondition(Fu)

model.bcs.T.top = FluxBoundaryCondition(Fθ)
model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

fig, axs = subplots(ncols=4, sharey=true)

sca(axs[1]); xlabel(L"\Delta T")
sca(axs[2]); xlabel(L"T")
sca(axs[3]); xlabel(L"\overline{w^2}")
sca(axs[4]); xlabel(L"\epsilon")

#iterate!(model, Δt)

ModularKPP.update_state!(model)
makeplot!(axs, model)

#makeplot!(axs, model)
