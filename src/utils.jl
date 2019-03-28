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
    return quote
        struct Solution <: AbstractSolution{$(nfields), Field}
            $(solfields...)
        end

        struct BoundaryConditions <: FieldVector{$(nfields), FieldBoundaryConditions}
            $(bcfields...)
        end

        struct Functionary <: FieldVector{$(nfields), Function}
            $(opfields...)
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

#=
diff_N_fields = [ :( $(Symbol(:N, name))::Function    ) for name in names ]
diff_K_fields = [ :( $(Symbol(:K, name))::Function    ) for name in names ]
diffop_fields = [ :( $(name)::Tridiagonal{T, A}       ) for name in names ]

struct DiffusiveOperator{T, A} <: FieldVector{$(nfields), Tridiagonal{T, A}}
    $(diffop_fields...)
end

function DiffusiveOperator(solution)
    lhs_array = []
    for s in solution
        T = eltype(s)
        A = arraytype(s)
        N = length(s)
        lhs_s = Tridiagonal{T, A}(zeros(N-1), zeros(N), zeros(N-1))
        push!(lhs_array, lhs_s)
    end
    lhs = DiffusiveOperator(lhs_array...)
    return lhs
end
=#


@def add_standard_model_fields begin
  timestepper :: TS
  grid        :: G
  equation    :: E
  solution    :: Solution
  bcs         :: BoundaryConditions
  clock       :: Clock{T}
end
