using Pkg; Pkg.activate(".."); Pkg.instantiate()

using 
  OceanTurb.KPP,
  PyPlot, 
  PyPlotPlus

# First, we define the model.
model = Model(nz=200, Lz=100.0, K0=1)

# Next we use a simple initial condition and iterate forward
Lz = model.grid.Lz
z0, dz = -Lz/2, Lz/20
c0(z) = exp(-(z-z0)^2 / 2dz^2)

# Set c to the function c0(z) --- enabled by some OceanTurb.jl syntactic sugar
model.solution.T = c0

# Initialize plotting
ax, fig = subplots()
xlabel("T")
ylabel("z")
title("No flux boundary conditions")
cornerspines()

plot(model.solution.T.data, zdata(model.solution.T), ".", label="initial condition")

# Run
dt = 0.01
nt = 1000
ni = 10

for i = 1:ni
    iterate!(model, dt, nt)
    plot(model.solution.T.data, zdata(model.solution.T))
end

legend()
show()
