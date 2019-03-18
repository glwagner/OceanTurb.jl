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

#fieldtype(::AbstractSolution{N, T}) where {N, T} = T
#arraytype(::AbstractSolution{N, T}) where {N, T<:Field{L, A}} where {L, A} = A
#eltype(::AbstractSolution{N, T}) where {N, T<:Field{L, A}} where {L, A} = eltype(A)

macro specify_solution(T, names...)
    nfields = length(names)
    soln_fields = [ :( $(name)::$(T)                    ) for name in names ]
      bc_fields = [ :( $(name)::FieldBoundaryConditions ) for name in names ]

      esc(
        quote
            struct Solution <: AbstractSolution{$(nfields), $(T)}
                $(soln_fields...)
            end

            struct BoundaryConditions <: FieldVector{$(nfields), FieldBoundaryConditions}
                $(bc_fields...)
            end
        end
      )
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
  timestepper::TS
  grid::G
  equation::E
  solution::Solution
  bcs::BoundaryConditions
  clock::Clock{T}
end
