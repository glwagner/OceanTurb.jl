using Pkg; Pkg.activate(".."); Pkg.instantiate()

using OceanTurb, Printf

@use_pyplot_utils

model = Diffusion.Model(
    N = 100,
    L = 1.0,
    K = 0.01,
    stepper = :BackwardEuler)

c₀(z) = exp(-(z + model.grid.L/2)^2/(2*(model.grid.L/20)^2))
model.solution.c = c₀

fig, axs = subplots()
xlabel(L"c")
ylabel(L"z")

plot(model.solution.c, label=L"t=0")
gcf()

iterate!(model, Δt=0.01, Nt=100)

plot(model.solution.c, label=@sprintf("\$ t = %0.2f \$", time(model)))
removespines("top", "right")

legend()
gcf()
