using JLD2, Printf, OceanTurb

# For this example we use the same "simple flux" scenario used in kpp_examples.ipynb.

function simple_flux_model(; N=10, L=100, Tz=0.01, Fb=1e-8, Fu=0, parameters=KPP.Parameters())

    model = KPP.Model(N=N, L=L, parameters=parameters, stepper=:BackwardEuler)

    # Initial condition
    T₀(z) = 20 + Tz*z

    # Set T to the function T0(z)
    model.solution.T = T₀

    # Fluxes
    Qθ = Qb / (model.constants.α * model.constants.g)
    model.bcs.U.top = FluxBoundaryCondition(Qu)
    model.bcs.T.top = FluxBoundaryCondition(Qθ)
    model.bcs.T.bottom = GradientBoundaryCondition(Tz)

    return model
end

## # Saving data
#
# To facilitate data, we define a few functions that save useful things both on model 
# initialization (init_data), and while the model is running (save_data).
# 
# We use JLD2.jl.

function save_data(path, model)
    iteration = iter(model)
    jldopen(path, "a+") do file
        file["timeseries/t/$iteration"] = time(model)
        file["timeseries/U/$iteration"] = collect(data(model.solution.U))
        file["timeseries/T/$iteration"] = collect(data(model.solution.T))
    end
end

function init_data(path, model)
    Fu = OceanTurb.getbc(model, model.bcs.U.top)
    Fb = OceanTurb.getbc(model, model.bcs.T.top) * model.constants.α * model.constants.g
    Tz = OceanTurb.getbc(model, model.bcs.T.bottom)
    
    jldopen(path, "a+") do file
        for gridfield in (:N, :L)
            file["grid/$gridfield"] = getproperty(model.grid, gridfield)
        end

        file["bcs/Fu"] = Fu
        file["bcs/Fb"] = Fb
        file["bcs/Tz"] = Tz
    end

    return nothing
end

# Next, we generate a model with some benign parameters and run, saving data every `dout` seconds.

dt = 10 # time step in seconds
dout = 10*minute # interval between output
tfinal = 4*hour # final time

ntot = Int(tfinal/dt)

nint = Int(dout/dt)
nout = Int(ntot/nint)

filepath = joinpath(".", "test_free_convection.jld2")

model = simple_flux_model(
     N = 100,
     L = 100,
    Tz = 4.08e-4,
    Fb = 5e-8,
    Fu = 0
)

isfile(filepath) && rm(filepath)
init_data(filepath, model)

for i = 1:nout
    @printf "\nrunning from iter %d" iter(model)
    @time iterate!(model, dt, nint)

    @printf "    ...saving at iter %d" iter(model)
    @time save_data(filepath, model)
end

# # Opening and analyzing data
#
# Running the model created a jld file, which we can open and look at:

file = jldopen(filepath, "r")

@show file

N = file["grid/N"]
L = file["grid/L"]
@show iters = parse.(Int, keys(file["timeseries/t"]))

close(file)

# # Plotting data
#
# Now we can plot the data (and otherwise analyze it if we wish).

grid = UniformGrid(N, L)

function get_t_and_T(path, i)
    file = jldopen(path, "r")
    t = file["timeseries/t/$i"]
    T = file["timeseries/T/$i"]
    close(file)
    t, T
end

fig, axs = subplots()
cornerspines()
xlabel(L"T")
ylabel(L"z")

for (ii, iter) in enumerate(iters)
    t, T = get_t_and_T(filepath, iter)
    plotlabel = ii % 4 == 0 ? @sprintf("\$ t = %.2f \$ hours", t/hour) : ""
    plot(T, grid.zc, label=plotlabel)
end

legend()
