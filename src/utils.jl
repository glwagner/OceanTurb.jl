import Base: eltype

# Gregorian calendar-ic globals
const second = 1.0
const minute = 60second
const hour = 60minute
const day = 24hour
const year = 365day
const stellaryear = 23hour + 56minute + 4.098903691
const Ω = 2π/stellaryear

@inline nan2inf(a::T) where T = ifelse(isnan(a), T(Inf), a)
@inline inf2zero(a::T) where T = ifelse(isinf(a), zero(T), a)

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

"""
    run_until!(model, dt, tfinal)

Run `model` until `tfinal` with time-step `dt`.
"""
function run_until!(model, dt, tfinal)
    nt = floor(Int, (tfinal - time(model))/dt)
    iterate!(model, dt, nt)

    last_dt = tfinal - time(model)
    iterate!(model, last_dt)

    return nothing
end

"""
    diffusive_flux!(flux, fieldname, model)

Calculates the diffusive flux associated with `fieldname` in `model.solution.fieldname`
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

"""
    prettytime(t)

Convert a floating point value `t` representing an amount of time in seconds to a more
human-friendly formatted string with three decimal places. Depending on the value of `t`
the string will be formatted to show `t` in nanoseconds (ns), microseconds (μs),
milliseconds (ms), seconds (s), minutes (min), hours (hr), or days (day).
"""
function prettytime(t)
    # Modified from: https://github.com/JuliaCI/BenchmarkTools.jl/blob/master/src/trials.jl
    iszero(t) && return "0.000 s"

    if t < 1e-6
        value, units = t * 1e9, "ns"
    elseif t < 1e-3
        value, units = t * 1e6, "μs"
    elseif t < 1
        value, units = t * 1e3, "ms"
    elseif t < minute
        value, units = t, "s"
    elseif t < hour
        value, units = t / minute, "min"
    elseif t < day
        value, units = t / hour, "hr"
    else
        value, units = t / day, "day"
    end

    return @sprintf("%.3f", value) * " " * units
end
