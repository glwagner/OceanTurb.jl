using OceanTurb, Printf 

using OceanTurb.ModularKPP: DiagnosticPlumeModel

@use_pyplot_utils

usecmbright()

function makeplot!(axs, model)
    TKEMassFlux.update_state!(model)

    ΔT = model.state.plume.T - model.solution.T

    ϵ = CellField(model.grid)
    for i in eachindex(ϵ)
        @inbounds ϵ[i] = TKEMassFlux.tracer_entrainment(model, i)
    end

    for ax in axs
        sca(ax); cla()
    end

    sca(axs[1]); xlabel(L"\breve T - T"); ylabel(L"z \, \mathrm{(m)}")
    sca(axs[2]); xlabel(L"T")
    sca(axs[3]); xlabel(L"\breve W^2")
    sca(axs[4]); xlabel(L"\epsilon")

    sca(axs[1])
    plot(ΔT)
    removespines("top", "right")

    sca(axs[2])
    plot(model.solution.T, label="Environment \$ T \$")
    plot(model.state.plume.T, "--", label="Plume \$ \\breve T \$")
    removespines("top", "left", "right")

    sca(axs[3])
    plot(model.state.plume.W²)
    removespines("top", "left", "right")

    sca(axs[4])
    plot(ϵ)
    removespines("top", "left", "right")

    axs[2].tick_params(left=false, labelleft=false)
    axs[3].tick_params(left=false, labelleft=false)
    axs[4].tick_params(left=false, labelleft=false)

    return nothing
end

Qb = 1e-7
N² = 1e-5
Δt = 1second

# The standard setup except with a plume model rather than a counter-gradient flux model.
model = TKEMassFlux.Model(         grid = UniformGrid(N=64, H=64), 
                          nonlocal_flux = TKEMassFlux.DiagnosticPlumeModel(CQ=4.0, Ca=0.1),
                          mixing_length = TKEMassFlux.SimpleMixingLength(Cᴸʷ=10.0),
                                stepper = :BackwardEuler
                         )

# Initial condition and fluxes
dTdz = N² / (model.constants.α * model.constants.g)
T₀(z) = 20 + dTdz*z
model.solution.T = T₀
OceanTurb.set!(model.state.plume.T, T₀)

Qθ = Qb / (model.constants.α * model.constants.g)

model.bcs.T.top = FluxBoundaryCondition(Qθ)
model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

close("all")

fig, axs = subplots(ncols=4, sharey=true)

run_until!(model, Δt, 1minute)
makeplot!(axs, model)
