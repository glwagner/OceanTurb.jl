# --
# Define tests
# --

function test_zeros(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  OceanTurb.@zeros T dims a2 b2
  eltype(a1) == T && a1 == a2 && b1 == b2
end
