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

function build_solution(fieldnames, fieldtypes=[Field for name in fieldnames]; name=Symbol(""))
    nfields = length(fieldnames)
    solfields = [ :( $(fieldnames[i]) :: $(fieldtypes[i]) ) for i = 1:nfields ]
    bcfields =  [ :( $(fieldnames[i]) :: FieldBoundaryConditions ) for i = 1:nfields ]
    opfields =  [ :( $(fieldnames[i]) :: Function ) for i = 1:nfields ]
    accfields = [ :( $(fieldnames[i]) :: T ) for i = 1:nfields ]
    lhsfields = [ :( $(fieldnames[i]) :: Tridiagonal{T, A} ) for i = 1:nfields ]

    sol_name = Symbol(name, :Solution)
    bc_name = Symbol(name, :BoundaryConditions)
    acc_signature = Expr(:curly, Symbol(name, :Accessory), :T)
    lhs_signature = Expr(:curly, Symbol(name, :LeftHandSide), :T, :A)

    return quote
        import StaticArrays: FieldVector
        import LinearAlgebra: Tridiagonal
        import OceanTurb: build_lhs

        struct $sol_name <: AbstractSolution{$(nfields), Field}
            $(solfields...)
        end

        struct $bc_name <: FieldVector{$(nfields), FieldBoundaryConditions}
            $(bcfields...)
        end

        struct $acc_signature <: FieldVector{$(nfields), T}
            $(accfields...)
        end

        struct $lhs_signature <: FieldVector{$(nfields), Tridiagonal{T, A}}
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

macro solution(fieldnames...)
    esc(build_solution(fieldnames))
end

macro named_solution(name, fieldnames...)
    esc(build_solution(fieldnames, name=name))
end

macro typed_solution(T, fieldnames...)
    fieldtypes = [ T for name in fieldnames ]
    esc(build_solution(fieldnames, fieldtypes))
end

macro pair_typed_solution(paired_specs...)
    fieldnames = Symbol[]
    fieldtypes = Symbol[]
    for (i, spec) in enumerate(paired_specs)
        isodd(i) && push!(fieldtypes, spec)
        iseven(i) && push!(fieldnames, spec)
    end
    esc(build_solution(fieldnames, fieldtypes))
end

@def add_standard_model_fields begin
    clock       :: Clock{T}
    grid        :: G
    timestepper :: TS
    solution    :: Solution
    bcs         :: BoundaryConditions
end

@def add_clock_grid_timestepper begin
    clock       :: Clock{T}
    grid        :: G
    timestepper :: TS
end

function run_until!(model, dt, tfinal)
    nt = floor(Int, tfinal/dt)
    iterate!(model, dt, nt)

    last_dt = time(model) - tfinal
    iterate!(model, last_dt)
    return nothing
end
