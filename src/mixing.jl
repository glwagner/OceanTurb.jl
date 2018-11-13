function diffuse!(u, κ)
  @. u *= (1-2κ)
  u[1] = u[2] # no flux bottom boundary condition
  @views @. u[2:end-1] = κ*(u[1:end-2] + u[2:end])
  nothing
end

function bulkmixing!(profile, Ri) 
  nothing
end
