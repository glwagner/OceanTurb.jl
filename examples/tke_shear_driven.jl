using OceanTurb

@use_pyplot_utils # add utilities for plotting OceanTurb Fields

constants = Constants(f=1e-4)

     N = 128        # Model resolution
     L = 128        # Vertical extent of the model domain
    Qᵘ = -1e-4      # Surface buoyancy flux (positive implies cooling)
    N² = 1e-5        # Interior/initial temperature gradient
    Δt = 1minute   # Time step size
tfinal = 1hour      # Final time

dTdz = constants.α * constants.g * N²

# Build the model with a Backward Euler timestepper
model = TKEMassFlux.Model(N=N, L=L, stepper=:BackwardEuler, constants=constants,
                          mixing_length=TKEMassFlux.SimpleMixingLength(CLz=0.4, CLb=0.1, CLΔ=1.0),
                          tke_equation=TKEMassFlux.TKEParameters(CDe=0.0305, CK_T=2.0, CK_e=0.01, CK_U=0.8))

# Set initial condition
T₀(z) = 20 + dTdz * z
model.solution.T = T₀

# Set boundary conditions
model.bcs.U.top = FluxBoundaryCondition(Qᵘ)
model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

# Run the model
run_until!(model, Δt, tfinal)

ℓ = FaceField(model.grid)
K = FaceField(model.grid)
for i in eachindex(ℓ)
    ℓ[i] = TKEMassFlux.mixing_length_face(model, i)
    K[i] = TKEMassFlux.KU(model, i)
end

fig, axs = subplots(ncols=5, figsize=(16, 6))

sca(axs[1])
plot(model.solution.T)
removespines("top", "right")
xlabel("Temperature (\$ {}^\\circ \\mathrm{C} \$)")
ylabel(L"z \, \mathrm{(m)}")

sca(axs[2])
plot(model.solution.e)
removespines("top", "right", "left")
xlabel(L"e")

sca(axs[3])
plot(ℓ)
removespines("top", "right", "left")
xlabel(L"\ell")

sca(axs[4])
plot(K)
removespines("top", "right", "left")
xlabel(L"K_U")

sca(axs[5])
plot(model.solution.U, label=L"U")
plot(model.solution.V, label=L"V")
removespines("top", "right", "left")
legend()
xlabel(L"U")

for i = 2:length(axs)
    axs[i].tick_params(left=false, labelleft=false)
end

#gcf()
#savefig("figs/kpp_free_convection.png", dpi=480)
