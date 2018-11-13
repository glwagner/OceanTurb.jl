using
  PriceWellerPinkel,
  Test

function testzeros(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  @zeros T dims a2 b2
  eltype(a1) == T && a1 == a2 && b1 == b2
end

function testexample(H=200, nz=200)
  examplemodel = loadexample(H, nz)
  nz == examplemodel.profile.nz
end

function unstablemodel(H=100, z1=20, z2=40, nz=100)
  examplemodel = loadexample(H, nz)
  p = examplemodel.profile
  for i = 1:nz
    if p.z[i] < z2
      p.ρ[i] = 2ρ₀
    elseif p.z[i] < z1
      p.ρ[i] = 0.5*ρ₀
    else
      p.ρ[i] = ρ₀
    end
  end
  nothing
end

@test testzeros(Float64)
@test testzeros(Float32)
@test testexample()
