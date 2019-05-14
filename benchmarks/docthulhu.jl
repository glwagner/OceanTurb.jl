using Cthulhu, OceanTurb
Δt = 1
N = 128
m = KPP.Model(N=N, stepper=:BackwardEuler)
#@descend OceanTurb.calc_diffusive_lhs!(Δt, m.timestepper.lhs, m.timestepper.eqn.K, m.solution, m)
@descend OceanTurb.calc_implicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.solution, m)
#OceanTurb.calc_diffusive_lhs!(Δt, m.timestepper.lhs, m.timestepper.eqn.K, m.solution, m)

