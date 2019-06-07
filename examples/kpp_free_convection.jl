using OceanTurb

@use_pyplot_utils

     N = 128
     L = 128
    Fb = 1e-7
  dTdz = 1e-3
    Δt = 10minute
tfinal = 8hour

model = KPP.Model(N=N, L=L, stepper=:BackwardEuler)
model.solution.T = T₀

model.bcs.T.top = FluxBoundaryCondition(Fb / (model.constants.α * model.constants.g))
model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

run_until!(model, Δt, tfinal)

plot(model.solution.T)
removespines("top", "right")
xlabel("Temperature (\$ {}^\\circ \\mathrm{C} \$)")
ylabel(L"z \, \mathrm{(m)}")
gcf()

savefig("kpp_free_convection.png", dpi=480)
