using Pkg; Pkg.activate(".."); Pkg.instantiate()

using OceanTurb, Printf

@use_pyplot_utils

# Instantiate a diffusion model with 100 grid points, height 1.0,
# diffusivity 0.01, and a backward Euler timestepper.
model = Diffusion.Model(
    N = 100,
    L = 1.0,
    K = 0.01,
    stepper = :BackwardEuler)

# Set the initial condition to a Gaussian cenetered at z=-0.5.
c₀(z) = exp(-(z + 0.5)^2 / 0.005)
model.solution.c = c₀

# Iterate the model forward
iterate!(model, Δt=0.01, Nt=100)

# Plot some results
fig, axs = subplots()
xlabel(L"c")
ylabel(L"z")

plot(c₀.(model.grid.zc), label=L"t=0")
gcf()

plot(model.solution.c, label=@sprintf("\$ t = %0.2f \$", time(model)))
removespines("top", "right")

legend()
gcf()
