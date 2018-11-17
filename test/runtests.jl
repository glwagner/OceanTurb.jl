using
  PriceWellerPinkel,
  Test

const year = PriceWellerPinkel.year

function testzeros(T=Float64, dims=(13, 45))
    a1, b1 = zeros(T, dims), zeros(T, dims)
    @zeros T dims a2 b2
    eltype(a1) == T && a1 == a2 && b1 == b2
end

function testexample(H=200, nz=200)
    examplemodel = loadexample(H, nz)
    nz == examplemodel.profile.nz
end

function testconstants()
    constants = Constants(f=π)
    isapprox(constants.f, π)
end


function testforcinginterpolant(fld=:shortwave)
  kwargs = Dict(:tdata=>[0, year], fld=>[0, 1.0])
  forcing = Forcing(; kwargs...)
  finterp = ForcingInterpolant(forcing)
  isapprox(getfield(finterp, fld)(0.6year), 0.6)
end

function testconvect1()
    H, nz = 1, 4
    model = Model(H=H, nz=nz)
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
    model = Model(H=H, nz=nz)
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

constants = Constants(f=2π, ρ₀=1.0)
model = Model(constants=constants, H=1, nz=4)

function teststepuv1()
    model.profile.U .= 1.0
    #            imix  dt  τˣ τʸ
    step_U!(model, 2, 1/8, 0, 0)  
    updatevars!(model)
    # u^2 + v^2 = 1 => u = 1 / sqrt(2)
    uvanswer = ones(model.profile.nz) / sqrt(2)
    isapprox(model.profile.u, uvanswer) && isapprox(model.profile.v, -uvanswer)
end

function teststepuv2(; τˣ=1.0, imix=3, nt=10)
    model.profile.U .= 0.0
    dt = 1e-6
    tf = nt*dt
    for i = 1:nt
        step_U!(model, imix, dt, τˣ, 0)  
    end
    Uanswer = 0*model.profile.U
    h = model.profile.z[imix]
    f = model.constants.f
    ρ₀ = model.constants.ρ₀
    @views @. Uanswer[imix:end] = τˣ/(ρ₀*h) * (1 - exp(-im*f*tf)) # distributed at top of domain.
    @show model.profile.U, Uanswer
    isapprox(model.profile.U, Uanswer, rtol=1e-9)
end



@testset "Utilities" begin
    @test testzeros(Float64)
    @test testzeros(Float32)
end

@testset "Forcing" begin
    for fld in fieldnames(ForcingInterpolant)
        @test testforcinginterpolant(fld)
    end
end

@testset "Model" begin
    @test testconstants()
    @test testexample()
    @test testconvect1()
    @test testconvect2()
    @test teststepuv1()
    @test teststepuv2()
end
