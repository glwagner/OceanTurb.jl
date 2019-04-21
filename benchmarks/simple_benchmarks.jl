using Pkg; Pkg.activate("..")

using OceanTurb, Printf, StaticArrays

Ns = (100, 1000,)
steppers = (:ForwardEuler,)
nts = (1000,)
nups = (100,)

struct FakeSolution <: FieldVector{1, Array{AbstractFloat, 1}}
    c :: Array{Float64, 1}
end

struct AnotherFakeSolution <: FieldVector{1, Array{AbstractFloat, 1}}
    c :: Array{Float64, 1}
    d :: Array{Float32, 1}
end




# Measure memory allocation
#=
@printf "\nTesting model instantiation...\n"
for stepper in steppers
    for N in Ns
        @printf "stepper: % 16s, N: % 6d" stepper N
        @time model = Diffusion.Model(N=N, stepper=stepper)
    end
end

function manyupdates!(m, n)
    for i = 1:n
        OceanTurb.update!(m)
    end
    return nothing
end

@printf "\nTesting updates...\n"
for N in Ns
    for stepper in steppers
        model = Diffusion.Model(N=N, stepper=stepper)
        manyupdates!(model, 1)
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            @time manyupdates!(model, nt)
        end
    end
end

function manyunpacks!(m, n)
    for j in eachindex(m.solution)
        for i = 1:n
            @inbounds ϕ, rhsϕ, Rϕ, Kϕ, bcsϕ = OceanTurb.unpack(m, j)
        end
    end
    return nothing
end

@printf "\nTesting unpack...\n"
for N in Ns
    for stepper in steppers
        model = Diffusion.Model(N=N, stepper=stepper)
        manyunpacks!(model, 1)
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            @time manyunpacks!(model, nt)
        end
    end
end
=#

function manyrhs!(m, n)
    for i = 1:n
        #OceanTurb.calc_explicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m)
        for j in eachindex(m.solution)
            ϕ = m.solution[j]
            #ϕ, rhsϕ, Rϕ, Kϕ, bcsϕ = OceanTurb.unpack(m, j)
            for i in eachindex(ϕ)
                @inbounds ϕ.data[i] = 0.0
                #rhsϕ[i] = OceanTurb.explicit_rhs_kernel(ϕ, Kϕ, Rϕ, m, i)
                #m.timestepper.rhs[j][i] = OceanTurb.explicit_rhs_kernel(ϕ, Kϕ, Rϕ, m, i)
            end
        end
    end
    return nothing
end

function manyfakerhs!(solution, n)
    for i = 1:n
        #OceanTurb.calc_explicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m)
        for j in eachindex(solution)
            ϕ = solution[j]
            for i in eachindex(ϕ)
                @inbounds ϕ[i] = 0.0
            end
        end
    end
    return nothing
end


@printf "\nTesting calc rhs...\n"
for N in Ns
    for stepper in (:ForwardEuler,)
        model = Diffusion.Model(N=N, stepper=stepper)
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            manyrhs!(model, nt)
            @time manyrhs!(model, nt)
        end
    end
end

@printf "\nTesting a fake calc rhs...\n"
for N in Ns
    for stepper in (:ForwardEuler,)
        c = rand(N)
        solution = FakeSolution(c)
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            manyfakerhs!(solution, nt)
            @time manyfakerhs!(solution, nt)
        end
    end
end

@printf "\nTesting another fake calc rhs...\n"
for N in Ns
    for stepper in (:ForwardEuler,)
        c = rand(N)
        d = rand(N)
        solution = AnotherFakeSolution(c, d)
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            manyfakerhs!(solution, nt)
            @time manyfakerhs!(solution, nt)
        end
    end
end





#=
@printf "\nTesting timestepping...\n"
for N in Ns
    for stepper in steppers
        model = Diffusion.Model(N=N, stepper=stepper)
        iterate!(model, 0.1)

        # Measure memory allocation
        for nt in nts
            @printf "stepper: % 16s, N: % 6d, nt: % 6d" stepper N nt
            @time iterate!(model, 0.1, nt)
        end
    end
end
=#
