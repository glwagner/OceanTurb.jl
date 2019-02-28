# --
# Define tests
# --

function test_diffusion_basic()
  model = Diffusion.Model(nz=4, Lz=2.0, κ=0.1)
  model.parameters.κ == 0.1
end

function test_diffusion_set_c()
  model = Diffusion.Model(nz=4, Lz=2.0, κ=0.1)
  c0 = [1.0, 2.0, 3.0, 4.0]
  Diffusion.set!(model, c=c0)
  model.solution.data[1:model.grid.nz] == c0
end

function test_diffusion()
  model = Diffusion.Model(nz=100, Lz=1.0, κ=0.1)
  model.parameters.κ == 0.1
end
