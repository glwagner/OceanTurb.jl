using OceanTurb

using OceanTurb.TKEMassFlux: TKEParameters

@use_pyplot_utils # add utilities for plotting OceanTurb Fields

function makeplot!(fig, axs, model)

    fig.suptitle("\$ t = $(prettytime(model.clock.time)) \$")

    ℓ = CellField(model.grid)
    P = CellField(model.grid)
    B = CellField(model.grid)
    ϵ = CellField(model.grid)
    τ = CellField(model.grid)

    for i in eachindex(ℓ)
        ℓ[i] = TKEMassFlux.mixing_length(model, i)
        P[i] = TKEMassFlux.production(model, i)
        B[i] = TKEMassFlux.buoyancy_flux(model, i)
        ϵ[i] = TKEMassFlux.dissipation(model, i)
        τ[i] = isfinite(TKEMassFlux.tke_time_scale(model, i)) ? TKEMassFlux.tke_time_scale(model, i) : 0
    end

    sca(axs[1])
    cla()
    plot(model.solution.T, "s-")
    removespines("top", "right")
    xlabel("Temperature (\$ {}^\\circ \\mathrm{C} \$)")
    ylabel(L"z \, \mathrm{(m)}")
    xlim(model.solution.T[128-20], model.solution.T[128])
    grid()

    sca(axs[2])
    cla()
    plot(model.solution.e, "s-") # / -model.bcs.U.top.condition, "s")
    removespines("top", "right", "left")
    xlabel(L"e")
    grid()

    sca(axs[3])
    cla()
    plot(P, "s-", label="production")
    plot(B, "o--", label="buoyancy flux")
    plot(-ϵ, "^:", label="dissipation")
    removespines("top", "right", "left")
    xlabel(L"\partial_t e")
    #legend()
    grid()
    maxlim = maximum(abs, xlim())
    xlim(-maxlim, maxlim)

    sca(axs[4])
    cla()
    plot(ℓ, "s-")
    removespines("top", "right", "left")
    xlabel(L"\ell")
    grid()

    sca(axs[5])
    cla()
    plot(τ, "s-")
    removespines("top", "right", "left")
    xlabel(L"\tau")
    grid()

    #=
    sca(axs[5])
    cla()
    plot(model.state.K, "s-")
    removespines("top", "right", "left")
    xlabel(L"K")
    grid()
    =#

    sca(axs[6])
    cla()
    plot(model.solution.U, label=L"U", "s-")
    #plot(model.solution.V, label=L"V", "s")
    removespines("top", "right", "left")
    legend()
    xlabel(L"U")
    grid()

    for i = 2:length(axs)
        axs[i].tick_params(left=false, labelleft=false)
    end

    ylim(-20, 1)
    pause(0.1)

    return ℓ
end

constants = Constants(f=0.0)

 N = 128        # Model resolution
 L = 128        # Vertical extent of the model domain
Qᵘ = -1e-4      # Surface buoyancy flux (positive implies cooling)
N² = 1e-5       # Interior/initial temperature gradient
Δt = 1second

dTdz = constants.α * constants.g * N²

# Build the model with a Backward Euler timestepper
model = TKEMassFlux.Model(           grid = UniformGrid(N=N, L=L), 
                                  stepper = :BackwardEuler,
                             tke_equation = TKEParameters(Cᴾʳₑ=1e-1),
                           tke_wall_model = nothing,
                                constants = constants)

# Set initial condition
T₀(z) = 20 + dTdz * z
model.solution.T = T₀

# Set boundary conditions
model.bcs.U.top = FluxBoundaryCondition(Qᵘ)
model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

# Run the model

OceanTurb.update!(model)

fig, axs = subplots(ncols=6, figsize=(12, 5), sharey=true)

#ℓ = makeplot!(fig, axs, model)

while true
    iterate!(model, Δt, 1)
    #OceanTurb.update!(model)
    makeplot!(fig, axs, model)
    println("i = $(model.clock.iter)")
    readline()
end
