function test_zeros(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  OceanTurb.@zeros T dims a2 b2
  eltype(a1) == T && a1 == a2 && b1 == b2
end

function test_run_until(stepper=:ForwardEuler)
    model = Diffusion.Model(N=3, L=1, K=1.0, stepper=stepper)

    dt, tfinal = 0.5, 1.3
    run_until!(model, dt, tfinal)

    return time(model) == tfinal
end

function test_diffusive_flux()
    model = Diffusion.Model(N=4, L=2, K=0.1)
    c0(z) = z
    model.solution.c = c0
    flux = OceanTurb.diffusive_flux(:c, model)
    return flux[2] == -0.1
end

@testset "Utils" begin
    @test test_zeros(Float64)
    @test test_zeros(Float32)
    @test test_diffusive_flux()
    for stepper in steppers
        @test test_run_until(stepper)
    end
end

