using
  OceanMixedLayerModels,
  OceanMixedLayerModels.PriceWellerPinkel,
  Test

# OMLM: OceanMixedLayerModels
# PWP:  Price-Weller-Pinkel
# MY:   Mellor-Yamada
# KPP:  Kappa Profile Parameterization

const OMLM = OceanMixedLayerModels
const PWP = OceanMixedLayerModels.PriceWellerPinkel # temporary fix
const year = OMLM.year

# --
# Test running section
# --

@testset "Basic" begin
  include("testbasic.jl")
  #@test testparams()
  #@test testexample()
  @test testzeros(Float64)
  @test testzeros(Float32)
end

@testset "Forcing" begin

  function testforcinginterpolant(fld=:shortwave)
    kwargs = Dict(:tdata=>[0, year], fld=>[0, 1.0])
    forcing = Forcing(; kwargs...)
    finterp = ForcingInterpolant(forcing)
    isapprox(getfield(finterp, fld)(0.6year), 0.6)
  end

  for fld in fieldnames(ForcingInterpolant)
    @test testforcinginterpolant(fld)
  end
end

#=
@testset "Mixing" begin
  include("testmixing.jl")
  @test testconvect1()
  @test testconvect2()
  @test testgradri()
  @test testbulkri()
  @test testbulkmix1()
  @test testbulkmix2()
  @test testgradmixing()
end

@testset "Time stepping" begin
  include("teststeps.jl")
  @test teststepuv1()
  @test teststepuv2()
  @test teststepforward()
end
=#
