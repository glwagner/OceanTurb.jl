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

function build_homogeneous_solution(names, typeprefix=:CellField; prefix=Symbol(""))
    nfields = length(names)
    fieldtype = Expr(:curly, typeprefix, :A, :G, :T)
    exprs = [ :( $(names[i]) :: $fieldtype ) for i = 1:nfields ]

    fullname = Symbol(prefix, :Solution)
    signature = Expr(:curly, fullname, :A, :G, :T)

    return quote
        #import StaticArrays: FieldVector

        struct $signature <: AbstractSolution{$(nfields), $fieldtype}
            $(exprs...)
        end
    end
end

function build_accessories(names; prefix=Symbol(""))
    nfields = length(names)
    bcsfields = [ :( $(names[i]) :: FieldBoundaryConditions ) for i = 1:nfields ]
    accfields = [ :( $(names[i]) :: T ) for i = 1:nfields ]
    lhsfields = [ :( $(names[i]) :: Tridiagonal{T, A} ) for i = 1:nfields ]

    bcs_signature = Symbol(prefix, :BoundaryConditions)
    acc_signature = Expr(:curly, Symbol(prefix, :Accessory), :T)
    lhs_signature = Expr(:curly, Symbol(prefix, :LeftHandSide), :T, :A)

    return quote
        import StaticArrays: FieldVector
        import LinearAlgebra: Tridiagonal
        import OceanTurb: build_lhs

        struct $bcs_signature <: FieldVector{$(nfields), FieldBoundaryConditions}
            $(bcsfields...)
        end

        struct $acc_signature <: FieldVector{$(nfields), T}
            $(accfields...)
        end

        struct $lhs_signature <: FieldVector{$(nfields), Tridiagonal{T, A}}
            $(lhsfields...)
            function LeftHandSide(solution::AbstractSolution{1, AbstractField})
                lhs = build_lhs(solution)[1]
                new{eltype(solution[1]), arraytype(solution[1])}(lhs)
            end
        end

        function LeftHandSide(solution::AbstractSolution{N, AbstractField}) where N
            lhs = build_lhs(solution)
            LeftHandSide{eltype(solution[1]), arraytype(solution[1])}(lhs...)
        end
    end
end

macro solution(names...)
    esc(quote
        $(build_homogeneous_solution(names))
        $(build_accessories(names))
    end
    )
end

macro prefixed_solution(prefix, names...)
    esc(quote
        $(build_homogeneous_solution(names; prefix=prefix))
        $(build_accessories(names; prefix=prefix))
    end
    )
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
    nt = floor(Int, (tfinal - time(model))/dt)
    iterate!(model, dt, nt)

    last_dt = tfinal - time(model)
    iterate!(model, last_dt)

    return nothing
end

"""
    diffusive_flux!(flux, fieldname, model)

Calculates the diffusive flux associated with `fieldname` in `model.solution.fieldname`i
and stores the result in the `FaceField` `flux`.
"""
function diffusive_flux!(flux, fieldname, model)
    field = getproperty(model.solution, fieldname)
    K = getproperty(model.timestepper.eqn.K, fieldname)

    for i in interiorindices(flux)
        @inbounds flux[i] = - K(model, i) * ∂z(field, i)
    end

    return nothing
end

"""
    diffusive_flux(fieldname, model)

Returns `flux::FaceField` with the diffusive flux associated with `fieldname` 
in `model.solution.fieldname`.
"""
function diffusive_flux(fieldname, model)
    flux = FaceField(model.grid)
    diffusive_flux!(flux, fieldname, model)
    return flux
end
