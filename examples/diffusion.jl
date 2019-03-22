using Pkg; Pkg.activate(".."); Pkg.instantiate()

using OceanTurb

model = Diffusion.Model(nz=100, Lz=1.0, κ=0.01)
stepper = Timestepper(:ForwardEuler, model.solution)
