using Pkg; Pkg.activate("..")

using OceanTurb, Printf, StaticArrays, BenchmarkTools

Ns = (32,)

function atomic_forward_euler_benchmark(m, Δt)
    @printf "% 24s:" "prepare"
    @btime OceanTurb.prepare!($(m.bcs), $(m.timestepper.eqn.K), $(m.solution), $(m.timestepper.eqn.update!), $(m.grid.N), $m)

    @printf "% 24s:" "calc explicit rhs"
    @btime OceanTurb.calc_explicit_rhs!($(m.timestepper.rhs), $(m.timestepper.eqn), $(m.solution), $m)

    @printf "% 24s:" "update solution"
    @btime OceanTurb.forward_euler_update!($(m.timestepper.rhs), $(m.solution), $Δt)

    @printf "% 24s:" "tick"
    @btime OceanTurb.tick!($(m.clock), $Δt)
    return nothing
end

function atomic_backward_euler_benchmark(m, Δt)
    @printf "% 24s:" "prepare"
    @btime OceanTurb.prepare!($(m.bcs), $(m.timestepper.eqn.K), $(m.solution), $(m.timestepper.eqn.update!), $(m.grid.N), $m)

    @printf "% 24s:" "calc implicit rhs"
    @btime OceanTurb.calc_implicit_rhs!($(m.timestepper.rhs), $(m.timestepper.eqn), $(m.solution), $m)

    @printf "% 24s:" "calc diffusive lhs"
    @btime OceanTurb.calc_diffusive_lhs!($Δt, $(m.timestepper.lhs), $(m.timestepper.eqn.K), $(m.solution), $m)
    
    @printf "% 24s:" "update solution"
    @btime OceanTurb.backward_euler_update!($(m.timestepper.rhs), $(m.timestepper.lhs), $(m.solution), $Δt)

    @printf "% 24s:" "tick"
    @btime OceanTurb.tick!($(m.clock), $Δt)
    return nothing
end


@printf "\nAtomic timestepping benchmarking of default configuration...\n"
Δt = 1.1
for N in Ns

    @printf "N: % 4d, forward Euler\n" N
    m = Diffusion.Model(N=N, stepper=:ForwardEuler)
    atomic_forward_euler_benchmark(m, Δt)

    @printf "N: % 4d, backward Euler\n" N
    m = Diffusion.Model(N=N, stepper=:BackwardEuler)
    atomic_backward_euler_benchmark(m, Δt)
end


@printf "\nAtomic timestepping benchmarking of configuration with one constant boundary condition...\n"
Δt = 1.1
for N in Ns

    @printf "N: % 4d, forward Euler\n" N
    m = Diffusion.Model(N=N, stepper=:ForwardEuler)
    m.bcs.c.top = BoundaryCondition(Gradient, 0)
    atomic_forward_euler_benchmark(m, Δt)

    @printf "N: % 4d, backward Euler\n" N
    m = Diffusion.Model(N=N, stepper=:BackwardEuler)
    m.bcs.c.top = BoundaryCondition(Gradient, 0)
    atomic_backward_euler_benchmark(m, Δt)
end


@printf "\nAtomic timestepping benchmarking of configuration with one function boundary condition...\n"
Δt = 1.1
bc(model) = time(model)
for N in Ns

    @printf "N: % 4d, forward Euler\n" N
    m = Diffusion.Model(N=N, stepper=:ForwardEuler)
    m.bcs.c.top = BoundaryCondition(Gradient, bc)
    atomic_forward_euler_benchmark(m, Δt)

    @printf "N: % 4d, backward Euler\n" N
    m = Diffusion.Model(N=N, stepper=:BackwardEuler)
    m.bcs.c.top = BoundaryCondition(Gradient, bc)
    atomic_backward_euler_benchmark(m, Δt)
end

