using Pkg; Pkg.activate("..")

using OceanTurb, Printf, StaticArrays, BenchmarkTools

Ns = (100,)
steppers = (:ForwardEuler,)
nts = (100,)

struct FakeSolution{T} <: FieldVector{1, AbstractArray}
    c :: Array{T, 1}
    function FakeSolution(a::Array{T, 1}) where T <: AbstractFloat
        new{T}(a)
    end
end

struct BadFakeSolution <: FieldVector{1, Array{AbstractFloat, 1}}
    c :: AbstractArray
end

struct AnotherFakeSolution{T} <: FieldVector{2, AbstractField}
    c :: T
    d :: T
end

struct NonVectorFakeSolution{A, G, T}
    c :: FaceField{A, G, T}
    d :: FaceField{A, G, T}
end

propertynames(::NonVectorFakeSolution) = (:c, :d)

function build_solution(name, fieldnames)
    nfields = length(fieldnames)
    fields = [ :( $(fieldnames[i]) :: A ) for i = 1:nfields ]

    sol_name = Symbol(name, :Solution)
    signature = Expr(:curly, sol_name, :A)

    if nfields == 1
        constructor = quote
            function $sol_name(a::A) where A <: AbstractArray
                new{A}(a)
            end
        end
    else
        constructor = ""
    end

    return quote
        import StaticArrays: FieldVector

        struct $signature <: AbstractSolution{$(nfields), AbstractArray}
            $(fields...)
            $constructor
        end
    end
end

macro solution(name, fieldnames...)
    esc(build_solution(name, fieldnames))
end

function manyrhs!(m, n)
    for i = 1:n
        OceanTurb.calc_explicit_rhs!(m.timestepper.rhs, m.timestepper.eqn, m)
    end
    return nothing
end

function manyfakerhs!(solution, n)
    for i = 1:n
        for j in eachindex(solution)
            @inbounds ϕ = data(solution[j])
            for i in eachindex(ϕ)
                @inbounds ϕ[i] = rand()
            end
        end
    end
    return nothing
end

function manyfakerhs!(solution::Array, n)
    for i = 1:n
        for j in eachindex(solution)
            @inbounds ϕ = solution[j]
            for i in eachindex(ϕ)
                @inbounds ϕ[i] = rand()
            end
        end
    end
    return nothing
end

function manyfakerhs_fieldloop!(solution, n)
    for i = 1:n
        for fld in propertynames(solution)
            ϕfld = getproperty(solution, fld)
            ϕ = data(ϕfld)
            for i in eachindex(ϕ)
                @inbounds ϕ[i] = rand()
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
            @btime manyrhs!($model, $nt)
        end
    end
end

@printf "\nTesting another fake calc rhs...\n"
for N in Ns
    fc(z) = rand()
    fd(z) = rand()
    g = UniformGrid(N, 1.0)
    c = CellField(fc, g)
    d = CellField(fd, g)
    solution = AnotherFakeSolution(c, d)
    for nt in nts
        @printf "N: % 6d, nt: % 6d" N nt
        manyfakerhs!(solution, nt)
        @btime manyfakerhs!($solution, $nt)
    end
end

@printf "\nTesting loop of flds for another fake calc rhs...\n"
for N in Ns
    fc(z) = rand()
    fd(z) = rand()
    g = UniformGrid(N, 1.0)
    c = FaceField(fc, g)
    d = FaceField(fd, g)
    solution = NonVectorFakeSolution(c, d)
    for nt in nts
        @printf "N: % 6d, nt: % 6d" N nt
        manyfakerhs_fieldloop!(solution, nt)
        @btime manyfakerhs_fieldloop!($solution, $nt)
    end
end

@printf "\nTesting another fake calc rhs of arrays...\n"
for N in Ns
    fc(z) = rand()
    fd(z) = rand()
    g = UniformGrid(N, 1.0)
    c = FaceField(fc, g)
    d = CellField(fd, g)
    solution = [c, d]
    for nt in nts
        @printf "N: % 6d, nt: % 6d" N nt
        manyfakerhs!(solution, nt)
        @btime manyfakerhs!($solution, $nt)
    end
end
