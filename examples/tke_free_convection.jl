using OceanTurb

using OceanTurb.TKEMassFlux: TKEParameters

@use_pyplot_utils # add utilities for plotting OceanTurb Fields

function makeplot!(fig, axs, model)

    fig.suptitle("\$ t = \$ $(prettytime(model.clock.time))")

    ℓ = CellField(model.grid)
    P = CellField(model.grid)
    B = CellField(model.grid)
    ϵ = CellField(model.grid)

    markerkwargs = Dict(:marker=>"s", :markersize=>1)

    for i in eachindex(ℓ)
        ℓ[i] = TKEMassFlux.mixing_length(model, i)
        P[i] = TKEMassFlux.production(model, i)
        B[i] = TKEMassFlux.buoyancy_flux(model, i)
        ϵ[i] = TKEMassFlux.dissipation(model, i)
    end

    sca(axs[1])
    cla()
    plot(model.solution.T; linestyle="-", markerkwargs...)
    removespines("top", "right")
    xlabel("Temperature (\$ {}^\\circ \\mathrm{C} \$)")
    ylabel(L"z \, \mathrm{(m)}")

    sca(axs[2])
    cla()
    plot(model.solution.e; linestyle="-", markerkwargs...)
    removespines("top", "right", "left")
    xlabel(L"e")

    sca(axs[3])
    cla()
    plot(P; linestyle="-", label="production", markerkwargs...)
    plot(B; linestyle="--", label="buoyancy flux", markerkwargs...)
    plot(-ϵ; linestyle=":", label="dissipation", markerkwargs...)
    removespines("top", "right", "left")
    xlabel(L"\partial_t e")
    legend(loc="lower right")
    maxlim = maximum(abs, xlim())
    xlim(-maxlim, maxlim)

    sca(axs[4])
    cla()
    plot(ℓ; linestyle="-", markerkwargs...)
    removespines("top", "right", "left")
    xlabel(L"\ell")

    sca(axs[5])
    cla()
    plot(model.state.K; linestyle="-", markerkwargs...)
    removespines("top", "right", "left")
    xlabel(L"K")

    for i = 2:length(axs)
        axs[i].tick_params(left=false, labelleft=false)
    end

    ylim(-64, 1)
    pause(0.1)

    return ℓ
end

constants = Constants(f=1e-4)

 N = 128        # Model resolution
 H = 128        # Vertical extent of the model domain
Qᵇ = 1e-7       # Surface buoyancy flux (positive implies cooling)
N² = 1e-5       # Interior/initial temperature gradient
Δt = 1minute

dTdz = N² / (constants.α * constants.g)
Qᶿ = Qᵇ / (constants.α * constants.g)

# Build the model with a Backward Euler timestepper
model = TKEMassFlux.Model(                  grid = UniformGrid(N=N, H=H), 
                                         stepper = :BackwardEuler,
                           convective_adjustment = TKEMassFlux.FluxProportionalConvectiveAdjustment(Cᴬ=10.0),
                                       constants = constants)

# Set initial condition
T₀(z) = 20 + dTdz * z
model.solution.T = T₀

# Set boundary conditions
model.bcs.T.top = FluxBoundaryCondition(Qᶿ)
model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

# Run the model
fig, axs = subplots(ncols=5, figsize=(16, 5), sharey=true)

for iplot = 1:12
    run_until!(model, Δt, iplot * 1hour)
    OceanTurb.update!(model)
    makeplot!(fig, axs, model)
end
