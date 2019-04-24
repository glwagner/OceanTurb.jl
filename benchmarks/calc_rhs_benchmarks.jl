using Pkg; Pkg.activate("..")

using OceanTurb, Printf, StaticArrays, BenchmarkTools

Ns = (100,)
steppers = (:ForwardEuler,)
nts = (100,)

function manyrhs!(m, n)
    for i = 1:n
        OceanTurb.calc_explicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m.solution, m)
    end
    return nothing
end

function bad_manyrhs!(m, n)
    s1 = m.solution
    s2 = m.timestepper.rhs
    for i = 1:n
        for j in eachindex(s1)
            ϕ1 = s1[j]
            ϕ2 = s2[j]
            for i in eachindex(ϕ1)
                @inbounds ϕ1[i] = 2*ϕ2[i]
            end
        end
    end
    return nothing
end

function good_functionbarrier_manyrhs!(m, n)
    s1 = m.solution
    s2 = m.timestepper.rhs
    for i = 1:n
        fakerhs!(s1, s2)
    end
    return nothing
end

function fakerhs!(s1, s2)
    for j in eachindex(s1)
        ϕ1 = s1[j]
        ϕ2 = s2[j]
        for i in eachindex(ϕ1)
            @inbounds ϕ1[i] = 2*ϕ2[i]
        end
    end
    return nothing
end

@printf "\nTesting OceanTurb calc rhs...\n"
for N in Ns
    for stepper in (:ForwardEuler,)
        model = Diffusion.Model(N=N, stepper=stepper)
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            manyrhs!(model, nt)
            @btime manyrhs!($model, $nt)
        end
    end
end

@printf "\nTesting a bad method to loop over rhs...\n"
for N in Ns
    for stepper in (:ForwardEuler,)
        model = Diffusion.Model(N=N, stepper=stepper)
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            bad_manyrhs!(model, nt)
            @btime bad_manyrhs!($model, $nt)
        end
    end
end


@printf "\nTesting a good function barrier method to loop over rhs...\n"
for N in Ns
    for stepper in (:ForwardEuler,)
        model = Diffusion.Model(N=N, stepper=stepper)
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            good_functionbarrier_manyrhs!(model, nt)
            @btime good_functionbarrier_manyrhs!($model, $nt)
        end
    end
end
