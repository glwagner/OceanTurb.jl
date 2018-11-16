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

function testconvect1()
  H, nz = 1, 4
  model = TestModel(H, nz)
  p = model.profile
  # Unstable profile.
  p.ρ[4] = 4
  p.ρ[3] = 2
  p.ρ[2] = 5
  p.ρ[1] = 5

  ρanswer = zeros(nz)
  ρanswer[4] = 3
  ρanswer[3] = 3
  ρanswer[2] = 5
  ρanswer[1] = 5

  convect!(p)
  p.ρ == ρanswer
end

function testconvect2()
  H, nz = 1, 4
  model = TestModel(H, nz)
  p = model.profile
  # Unstable profile.
  p.ρ[4] = 6
  p.ρ[3] = 2
  p.ρ[2] = 2
  p.ρ[1] = 2

  p.T .= 1:4 # 1+2+3+4 = 10
  Tanswer = 2.5*fill(1, (nz,))

  ρanswer = zeros(nz)
  ρanswer[4] = 3
  ρanswer[3] = 3
  ρanswer[2] = 3
  ρanswer[1] = 3

  imix = convect!(p)

  isapprox(p.ρ, ρanswer) && isapprox(Tanswer, p.T)
end



@test testzeros(Float64)
@test testzeros(Float32)
@test testexample()
@test testconvect1()
@test testconvect2()
