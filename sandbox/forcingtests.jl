#=
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
=#
