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
