params = Parameters(f=2π, ρ₀=1.0, g=1.0, gradRiᵐⁱˣ=0.5)
model = PriceWellerPinkelModel(params=params, ocean=Ocean(H=4, nz=4))
ocean = model.ocean

function testconvect1()
  # Unstable ocean.
  ocean.ρ[4] = 4
  ocean.ρ[3] = 2
  ocean.ρ[2] = 5
  ocean.ρ[1] = 5

  ρanswer = zeros(ocean.nz)
  ρanswer[4] = 3
  ρanswer[3] = 3
  ρanswer[2] = 5
  ρanswer[1] = 5

  PWP.convect!(ocean, ocean.nz)
  ocean.ρ == ρanswer
end

function testconvect2()
  # Unstable ocean.
  ocean.ρ[4] = 6
  ocean.ρ[3] = 2
  ocean.ρ[2] = 2
  ocean.ρ[1] = 2

  ocean.T .= 1:4 # 1+2+3+4 = 10
  Tanswer = 2.5*fill(1, (ocean.nz,))

  ρanswer = zeros(ocean.nz)
  ρanswer[4] = 3
  ρanswer[3] = 3
  ρanswer[2] = 3
  ρanswer[1] = 3

  PWP.convect!(ocean, ocean.nz)
  isapprox(ocean.ρ, ρanswer) && isapprox(ocean.T, Tanswer)
end

function testgradri()
  # Simple example:
  U = ocean.U 
  ρ = ocean.ρ 
  model.imix = 3

  U[4] = 2
  U[3] = 2 # imix=3
  U[2] = 0
  U[1] = 0

  ρ[4] = 1
  ρ[3] = 1
  ρ[2] = 2
  ρ[1] = 2

  PWP.gradRi!(model)
  # Ri[4] = Inf
  # Ri[3] = 1/4
  # Ri[2] = Inf
  # Ri[1] = Inf

  model.ocean.Ri == [Inf, Inf, 0.25, Inf]
end

function testbulkri()
  # Simple example:
  U = ocean.U 
  ρ = ocean.ρ 
  model.imix = 3

  U[4] = 2 # z=0
  U[3] = 2 # z=-1, imix=3, h=1=-z[imix]
  U[2] = 0 # z=-2
  U[1] = 0 # z=-3

  ρ[4] = 1
  ρ[3] = 1
  ρ[2] = 2
  ρ[1] = 2

  h = 2
  Δρ = -1
  ΔU = 2
  Riᵇ = -h*Δρ/ΔU^2 
  
  PWP.bulkRi(model) == Riᵇ
end

function testbulkmix1()
  # Simple example:
  U = model.ocean.U 
  ρ = model.ocean.ρ 
  model.imix = 3

  U[4] = 1 # z=0
  U[3] = 1 # z=-1, imix=3, h=1=-z[imix]
  U[2] = 0 # z=-2
  U[1] = 0 # z=-3

  ρ[4] = 1
  ρ[3] = 1
  ρ[2] = 2
  ρ[1] = 2

  Uᶠ = deepcopy(U)
  ρᶠ = deepcopy(ρ)
  PWP.bulkmix!(model) # should do nothing

  model.ocean.U == Uᶠ && model.ocean.ρ == ρᶠ
end

function testbulkmix2()
  # Simple example:
  U = ocean.U 
  ρ = ocean.ρ 
  model.imix = 3

  U[4] = 2 # z=0
  U[3] = 2 # z=-1, imix=3, h=1=-z[imix]
  U[2] = 0 # z=-2
  U[1] = 0 # z=-3

  ρ[4] = 1
  ρ[3] = 1
  ρ[2] = 2
  ρ[1] = 2

  Uᶠ = 0U
  Uᶠ[4] = 4/3
  Uᶠ[3] = 4/3
  Uᶠ[2] = 4/3
  Uᶠ[1] = 0

  ρᶠ = 0ρ
  ρᶠ[4] = 4/3
  ρᶠ[3] = 4/3
  ρᶠ[2] = 4/3
  ρᶠ[1] = 2

  PWP.bulkmix!(model)

  model.ocean.U == Uᶠ && model.ocean.ρ == ρᶠ
end

function testgradmixing()
  # Simple example:
  U = model.ocean.U 
  ρ = model.ocean.ρ 
  model.imix = 3

  U[4] = 2 # z=0
  U[3] = 2 # z=-1, imix=3, h=1=-z[imix]
  U[2] = 0 # z=-2
  U[1] = 0 # z=-3

  ρ[4] = 1
  ρ[3] = 1
  ρ[2] = 2
  ρ[1] = 2

  # δU = 0.5*(1-0.5)*2 = -0.5
  Uᶠ = 0U
  Uᶠ[4] = 2
  Uᶠ[3] = 1.5
  Uᶠ[2] = 0.5
  Uᶠ[1] = 0

  # δρ = 0.5*(1-0.5)*1 = 0.25
  ρᶠ = 0ρ
  ρᶠ[4] = 1
  ρᶠ[3] = 1.25
  ρᶠ[2] = 1.75
  ρᶠ[1] = 2

  PWP.gradientmix!(model)

  model.ocean.U == Uᶠ && model.ocean.ρ == ρᶠ
end
