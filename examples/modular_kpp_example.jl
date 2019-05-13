using Pkg; Pkg.activate(".")

using OceanTurb, OceanTurb.Plotting, PyPlot, Printf

modelsetup = (N=100, L=100, stepper=:BackwardEuler)

Fb = 2e-8
Fu = 0.0
Tz = 0.001
Δt = 10*minute
tf = 8*hour

cvmix = KPP.Model(; modelsetup...)

holtslag = ModularKPP.Model(; modelsetup...,
    diffusivity=ModularKPP.HoltslagDiffusivityParameters())

roms = ModularKPP.Model(; modelsetup...,
    mixingdepth=ModularKPP.ROMSMixingDepthParameters())

# Initial condition and fluxes
T₀(z) = 20 + Tz*z

models = (cvmix, holtslag, roms)

for model in models
    model.solution.T = T₀

    Fθ = Fb / (model.constants.α * model.constants.g)
    model.bcs.U.top = FluxBoundaryCondition(Fu)
    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    model.bcs.T.bottom = GradientBoundaryCondition(Tz)
end

fig, axs = subplots()

removespines("top", "right")
xlabel("Temperature")
ylabel(L"z \, \mathrm{(m)}")

for i = 1:5
    if i > 1
        for model in models
            run_until!(model, Δt, (i-1)*12hour)
        end
    end

    @printf("""
        t : %.1f hours

          mixing depths
          =============

         cvmix h : %.2f
        holtslag h : %.2f
            roms h : %.2f
    \n""", time(cvmix)/hour,
    cvmix.state.h,
    holtslag.state.h,
    roms.state.h,
    )

    if i == 1
        vlabel = "CVMix KPP"
        hlabel = "Holtslag \$K\$-profile and CVMix mixing depth model"
        rlabel = "CVMix \$K\$-profile and ROMS mixing depth model"
    else
        vlabel = ""
        hlabel = ""
        rlabel = ""
    end

    if i == 1
        tlabel = text(cvmix.solution.T[end], 0,
            @sprintf("\$ t = %.0f \$ hours", time(cvmix)/hour),
            verticalalignment="bottom", horizontalalignment="center", color=defaultcolors[i])
    else
        tlabel = text(maximum(cvmix.solution.T.data), -cvmix.state.h,
            @sprintf("\$ t = %.0f \$ hours", time(cvmix)/hour),
            verticalalignment="bottom", horizontalalignment="left", color=defaultcolors[i])
    end

    plot(cvmix.solution.T,  "-",  color=defaultcolors[i], label=vlabel)
    plot(holtslag.solution.T, "--", color=defaultcolors[i], label=hlabel)
    plot(roms.solution.T,     ":",  color=defaultcolors[i], label=rlabel)
end

legend(fontsize=10)
gcf()

#=
fig, axs = subplots()

removespines("top", "right")

plot(cvmix.solution.T)

for model in models
    run_until!(model, Δt, tf)
end

plot(cvmix.solution.T)
plot(holtslag.solution.T)
plot(roms.solution.T)
gcf()
=#

#=
fig, axs = subplots()

removespines("top", "right")

plot(cvmix.solution.T)

for model in models
    run_until!(model, Δt, tf)
end

plot(cvmix.solution.T)
plot(holtslag.solution.T)
plot(roms.solution.T)
gcf()
=#
