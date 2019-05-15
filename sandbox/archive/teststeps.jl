params = Parameters(f=2π, ρ₀=1.0, g=1.0, gradRiᵐⁱˣ=0.5) # simple parameters for testing
model = PriceWellerPinkelModel(params=params, ocean=Ocean(H=1, nz=4)) # model with zero forcing

function teststepuv1()
  model.ocean.U .= 1.0
  model.imix = 2
  #                   dt  τˣ τʸ
  PWP.step_U!(model, 1/8, 0, 0)  
  updatevars!(model)
  # u^2 + v^2 = 1 => u = 1 / sqrt(2)
  uvanswer = ones(model.ocean.nz) / sqrt(2)
  isapprox(model.ocean.u, uvanswer) && isapprox(model.ocean.v, -uvanswer)
end


"""
Solve:

Ut + ifU = τˣ
=> ∂t (e^{ift} U) = e^{ift}*Gz
if U(0) = 0 => e^{ift} U = (1/if) * (e^{ift}-1)*τˣ
            => U = -(i/f) * (1 - e^{ift})*τˣ
"""
function teststepuv2(; τˣ=1.0, imix=3, nt=100)
  model.ocean.U .= 0.0
  model.imix = imix
  h = mixedlayerdepth(model.ocean.zᴳ, imix)
  f, ρ₀ = model.params.f, model.params.ρ₀
  dt = 1e-9
  tf = nt*dt
  for i = 1:nt
      PWP.step_U!(model, dt, τˣ, 0)  
  end
  Uanswer = 0*model.ocean.U
  @views @. Uanswer[imix:end] = im*τˣ/(f*ρ₀*h) * (1 - exp(-im*f*tf)) # distributed at top of domain.
  isapprox(model.ocean.U, Uanswer, rtol=1e-6)
end

function teststepforward()
  model.ocean.U .= 1.0
  model.imix = 2
  #                       dt
  PWP.stepforward!(model, 1/8)
  updatevars!(model)
  # u^2 + v^2 = 1 => u = 1 / sqrt(2)
  uvanswer = ones(model.ocean.nz) / sqrt(2)
  isapprox(model.ocean.u, uvanswer) && isapprox(model.ocean.v, -uvanswer)
end
