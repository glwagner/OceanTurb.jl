using Pkg; Pkg.activate("..")

using OceanTurb, Printf, StaticArrays, BenchmarkTools

Ns = (100,)
steppers = (:ForwardEuler,)
nts = (100,)

@printf "\nTesting OceanTurb calc rhs...\n"
for N in Ns
    for stepper in (:ForwardEuler,)
        m = Diffusion.Model(N=N, stepper=stepper)
        #m.bcs.c.top = OceanTurb.ConstantGradientBoundaryCondition(1)
        m.bcs.c.bottom = OceanTurb.GradientBoundaryCondition(1)
        Δt = 0.1
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d\n" stepper N nt

            @btime OceanTurb.prepare!($(m.bcs), $(m.timestepper.eqn.K), $(m.solution), $(m.timestepper.eqn.update!), $(m.grid.N), $m)
            @btime OceanTurb.calc_explicit_rhs!($(m.timestepper.rhs), $(m.timestepper.eqn), $(m.solution), $m)
            @btime OceanTurb.forward_euler_update!($(m.timestepper.rhs), $(m.solution), $Δt)
            @btime OceanTurb.tick!($(m.clock), $Δt)

        end
    end
end
