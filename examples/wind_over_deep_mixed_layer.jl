# # Wind over deep mixed layer
#
# In this example, we blow an inertially-oscillating atmospheric wind
# over a deep, pre-mixed boundary layer. We compare the evolution of the
# boundary layer with the K-Profile parameterization and the TKEMassFlux model.

using Plots, Printf

using OceanTurb

using OceanTurb.TKEMassFlux: ModelBoundaryConditions # This function works for TKE and KPP

# ## Resolution and domain size
#
# We use a grid with 256 points and 1 meter grid spacing.

grid = UniformGrid(N=256, H=256)


# ## Coriolis
constants = Constants(f=1e-4) # s⁻¹

# ## The inertial wind
#
# We blow an inertially-rotating wind with wind stress 1e-4 m² s⁻².
# This corresponds roughly to 6 m s⁻¹ atmospheric wind speed.

τ = 1e-4 # m² s⁻²

## Inertially rotating means Qᵘ + i Qᵛ = e⁻ⁱᶠᵗ τ:
Qᵘ_inertial(t, f) =   τ * cos(f * t)
Qᵛ_inertial(t, f) = - τ * sin(f * t)

## Boundary condition functions
Qᵘ(model) = Qᵘ_inertial(model.clock.time, model.constants.f)
Qᵛ(model) = Qᵛ_inertial(model.clock.time, model.constants.f)

U_boundary_conditions = BoundaryConditions(top = FluxBoundaryCondition(Qᵘ))
V_boundary_conditions = BoundaryConditions(top = FluxBoundaryCondition(Qᵛ))

boundary_conditions = ModelBoundaryConditions(U = U_boundary_conditions,
                                              V = V_boundary_conditions)

# We are now ready to build both models.

tke_model = TKEMassFlux.Model(     grid = grid,
                                stepper = :BackwardEuler,
                                    bcs = boundary_conditions,
                              constants = constants)

kpp_model = ModularKPP.Model(     grid = grid,
                               stepper = :BackwardEuler,
                                   bcs = boundary_conditions,
                             constants = constants)

# ## Stratification
#
# Finally, we prescribe a linear stratification with buoyancy frequency

N² = 5e-6 # s⁻²

# The temperature gradient given the default thermal expansion coefficient `α`
# and gravitational acceleration `g` is

dTdz = N² / (constants.α * constants.g)

# ## Initial mixed layer depth
#
# We prescribe an initial mixed layer depth

h = 100 # m

# such that the temperature gradient is a step function,

T₀(z) = 20 + dTdz * min(z + h, 0)

# We then set initial conditions and the temperature boundary condition for both models,

for model in (tke_model, kpp_model)
    model.solution.T = T₀
    model.bcs.T.bottom = GradientBoundaryCondition(dTdz)
end

# Finally, a few shortcuts make plotting a little easier.
T_tke = tke_model.solution.T
U_tke = tke_model.solution.U
V_tke = tke_model.solution.V

T_kpp = kpp_model.solution.T
U_kpp = kpp_model.solution.U
V_kpp = kpp_model.solution.V

z = nodes(U_tke)

speed(U, V) = sqrt.(data(U).^2 .+ data(V).^2)

max_speed = 0.6

# ## Run the model
#
# We use a 1 minute time-step, which is quite small and safe for these models.

Δt = 1minute

# Animate!

anim = @animate for iplot = 1:201

    time_label = Plots.text(@sprintf("t = %.2f hours", tke_model.clock.time / hour))

    temperature_plot = plot(data(T_tke), z,
                            legend = :topleft,
                             label = "TKE",
                            xlabel = "Temperature (ᵒC)",
                            ylabel = "z (m)",
                             xlims = (19.7, 20.05),
                           )

    speed_plot = plot(speed(U_tke, V_tke), z,
                           legend = false,
                           xlabel = "Speed (m/s)",
                            xlims = (0, max_speed),
                      annotations = (0.4, -200, time_label),
                     )

    velocities_plot = plot(data(U_tke), z,
                                legend = :bottomright,
                                 label = "U (TKE)",
                                 color = :blue,
                                xlabel = "U, V (m/s)",
                                 xlims = (-max_speed, max_speed),
                          )

    plot!(velocities_plot, data(V_tke), z, label="V (TKE)", linecolor=:red)
    plot!(velocities_plot, data(U_kpp), z, label="U (KPP)", linecolor=:blue, linestyle=:dash)
    plot!(velocities_plot, data(V_kpp), z, label="V (KPP)", linecolor=:red, linestyle=:dash)

    plot!(temperature_plot, data(T_kpp), z, label = "KPP")
    plot!(speed_plot, speed(U_kpp, V_kpp), z)

    plot(temperature_plot, speed_plot, velocities_plot, layout = (1, 3), size=(800, 500))

    run_until!(tke_model, Δt, iplot * hour / 2)
    run_until!(kpp_model, Δt, iplot * hour / 2)
end

mp4(anim, "wind_over_deep_mixed_layer.mp4", fps = 15)
