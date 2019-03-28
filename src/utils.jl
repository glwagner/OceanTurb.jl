import Base: eltype

# Gregorian calendar-ic globals
const second = 1.0
const minute = 60second
const hour = 60minute
const day = 24hour
const year = 365day
const stellaryear = 23hour + 56minute + 4.098903691
const Ω = 2π/stellaryear

function pressenter()
    println("\nPress enter to continue.")
    chomp(readline())
end

macro zeros(T, dims, vars...)
    expr = Expr(:block)
    append!(expr.args, [:($(esc(var)) = zeros($(esc(T)), $(esc(dims))); ) for var in vars])
    expr
end

macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

function build_solution(names, fieldtypes)
    nfields = length(names)
    solfields = [ :( $(names[i]) :: $(fieldtypes[i]) ) for i = 1:nfields ]
    bcfields =  [ :( $(names[i]) :: FieldBoundaryConditions ) for i = 1:nfields ]
    opfields =  [ :( $(names[i]) :: Function ) for i = 1:nfields ]
    ancfields = [ :( $(names[i]) :: T ) for i = 1:nfields ]
    lhsfields = [ :( $(names[i]) :: Tridiagonal{T, A} ) for i = 1:nfields ]
    return quote
        import StaticArrays: FieldVector
        import LinearAlgebra: Tridiagonal
        import OceanTurb: build_lhs

        struct Solution <: AbstractSolution{$(nfields), Field}
            $(solfields...)
        end

        struct BoundaryConditions <: FieldVector{$(nfields), FieldBoundaryConditions}
            $(bcfields...)
        end

        struct SolutionLike{T} <: FieldVector{$(nfields), T}
            $(ancfields...)
        end

        struct LeftHandSide{T, A} <: FieldVector{$(nfields), Tridiagonal{T, A}}
            $(lhsfields...)
            function LeftHandSide(solution::AbstractSolution{1, Field})
                lhs = build_lhs(solution)[1]
                new{eltype(solution[1]), arraytype(solution[1])}(lhs)
            end
        end

        function LeftHandSide(solution::AbstractSolution{N, Field}) where N
            lhs = build_lhs(solution)
            LeftHandSide{eltype(solution[1]), arraytype(solution[1])}(lhs...)
        end
    end
end

macro specify_solution(T, names...)
    fieldtypes = [ T for name in names ]
    esc(build_solution(names, fieldtypes))
end

macro pair_specify_solution(paired_specs...)
    names = Symbol[]
    fieldtypes = Symbol[]
    for (i, spec) in enumerate(paired_specs)
        isodd(i) && push!(fieldtypes, spec)
        iseven(i) && push!(names, spec)
    end
    esc(build_solution(names, fieldtypes))
end

@def add_standard_model_fields begin
  clock       :: Clock{T}
  grid        :: G
  timestepper :: TS
  solution    :: Solution
  bcs         :: BoundaryConditions
end
