export
  second,
  minute,
  hour,
  day,
  year,
  stellaryear,
  Ω,

  pressenter,
  @zeros

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
