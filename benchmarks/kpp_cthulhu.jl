using Cthulhu, OceanTurb

Δt = 1
N = 128
m_be = KPP.Model(N=N, stepper=:BackwardEuler)
m_fe = KPP.Model(N=N, stepper=:ForwardEuler)

m_fe.bcs.T.top = BoundaryCondition(Flux, 1.0)
m_be.bcs.T.top = BoundaryCondition(Flux, 1.0)

#@descend OceanTurb.calc_diffusive_lhs!(Δt, m.timestepper.lhs, m.timestepper.eqn.K, m.solution, m)
#@descend OceanTurb.calc_explicit_rhs!(m_fe.timestepper.rhs, m_fe.timestepper.eqn, m_fe.solution, m_fe)

OceanTurb.update!(m_be.bcs, m_be.timestepper.eqn, m_be.solution, m_be)

#@descend OceanTurb.calc_implicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.solution, m)
