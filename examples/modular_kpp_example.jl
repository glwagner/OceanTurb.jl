using Pkg; Pkg.activate(".")

using OceanTurb, OceanTurb.Plotting, PyPlot, Printf

modelsetup = (N=100, L=100, stepper=:BackwardEuler)

Fb = 2e-8
Fu = 0.0
Tz = 0.001
Δt = 10*minute
tf = 8*hour

vanilla = KPP.Model(; modelsetup...)

holtslag = ModularKPP.Model(; modelsetup...,
    diffusivity=ModularKPP.HoltslagDiffusivityParameters())

roms = ModularKPP.Model(; modelsetup...,
    mixingdepth=ModularKPP.ROMSMixingDepthParameters())

# Initial condition and fluxes
T₀(z) = 20 + Tz*z

models = (vanilla, holtslag, roms)

for model in models
    model.solution.T = T₀

    Fθ = Fb / (model.constants.α * model.constants.g)
    model.bcs.U.top = FluxBoundaryCondition(Fu)
    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    model.bcs.T.bottom = GradientBoundaryCondition(Tz)
end

fig, axs = subplots(ncols=2, sharey=true)

removespines("top", "right", ax=axs[1])

for ax in axs[2:end]
    removespines("left", "top", "right", ax=ax)
    ax.tick_params(left=false)
end

#sca(axs[1])
#plot(vanilla.solution.T)

#uke = FaceField(roms.grid)
#ker = FaceField(roms.grid)

for i = 1:4
    for model in models
        iterate!(model, Δt, 100)
    end
    @printf("""
        i : %d
         vanilla h : %.2f
        holtslag h : %.2f
            roms h : %.2f
    \n""", i,
    vanilla.state.h,
    holtslag.state.h,
    roms.state.h,
    )

    for i in eachindex(uke)
        uke[i] = - KPP.unresolved_kinetic_energy(roms, i) / roms.grid.zf[i]
        ker[i] = ModularKPP.h_kernel(roms, i)
    end

    sca(axs[1])
    plot(vanilla.solution.T, "-")
    plot(holtslag.solution.T, "--")
    plot(roms.solution.T, ":")

    sca(axs[2])
    plot(roms.state.h_criterion)
    xlabel("h criterion")
end

gcf()

#=
fig, axs = subplots()

removespines("top", "right")

plot(vanilla.solution.T)

for model in models
    run_until!(model, Δt, tf)
end

plot(vanilla.solution.T)
plot(holtslag.solution.T)
plot(roms.solution.T)
gcf()
=#
