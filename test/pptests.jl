#
# Tests for the Diffusion module
#

function test_pp_basic()
    model = PacanowskiPhilander.Model(N=4, L=2)
    model.grid.N == 4
end

function test_pp_diffusion_cosine()
    parameters = PacanowskiPhilander.Parameters(Cκ₀=1.0, Cκ₁=0.0)
    model = PacanowskiPhilander.Model(N=100, L=π/2, parameters=parameters)
    z = model.grid.zc

    c_init(z) = cos(2z)
    c_ans(z, t) = exp(-4t) * c_init(z)

    model.solution.T = c_init

    dt = 1e-3
    iterate!(model, dt)

    # The error tolerance is a bit arbitrary.
    norm(c_ans.(z, time(model)) .- data(model.solution.T)) < model.grid.N*1e-6
end
