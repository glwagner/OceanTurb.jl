# 
# Tests for the Diffusion module
# 

function test_diffusion_basic()
  model = Diffusion.Model(nz=4, Lz=2.0, κ=0.1)
  model.parameters.κ == 0.1
end

function test_diffusion_set_c()
  model = Diffusion.Model(nz=4, Lz=2.0, κ=0.1)
  c0 = [1.0, 2.0, 3.0, 4.0]
  model.solution.c = c0
  model.solution.c.data[1:model.grid.nz] == c0
end
