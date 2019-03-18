#
# Tests for the Diffusion module
#

function test_diffusion_basic()
    model = Diffusion.Model(N=4, L=2, κ=0.1)
    model.parameters.κ == 0.1
end

function test_diffusion_set_c()
    model = Diffusion.Model(N=4, L=2, κ=0.1)
    c0 = 1:4
    model.solution.c = c0
    model.solution.c.data[1:model.grid.N] == c0
end

function test_diffusion_cosine()
    model = Diffusion.Model(N=100, L=π/2, κ=1)
    z = model.grid.zc

    c_init(z) = cos(2z)
    c_ans(z, t) = exp(-4t) * c_init(z)

    model.solution.c = c_init

    dt = 1e-3
    iterate!(model, dt)

    # The error tolerance is a bit arbitrary.
    norm(c_ans.(z, model.clock.time) .- model.solution.c.data) < model.grid.N*1e-6
end
