using Cthulhu, OceanTurb

Δt = 1
N = 128
m_fe = Diffusion.Model(N=N, stepper=:ForwardEuler)

@descend OceanTurb.calc_explicit_rhs!(m_fe.timestepper.rhs, m_fe.timestepper.eqn, m_fe.solution, m_fe)
