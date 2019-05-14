using Pkg; Pkg.activate(".."); Pkg.instantiate()

using OceanTurb, Printf

@use_pyplot_utils

usecmbright()

modelsetup = (N=100, L=100, stepper=:BackwardEuler)

Fb = 2e-8
Fu = 0.0
Tz = 0.001
Δt = 10*minute
tf = 8*hour

        cvmix = ModularKPP.Model(; modelsetup...)

     holtslag = ModularKPP.Model(; modelsetup..., diffusivity = ModularKPP.HoltslagDiffusivityParameters())

         roms = ModularKPP.Model(; modelsetup..., mixingdepth = ModularKPP.ROMSMixingDepthParameters())

holtslag_roms = ModularKPP.Model(; modelsetup..., diffusivity = ModularKPP.HoltslagDiffusivityParameters(),
                                                  mixingdepth = ModularKPP.ROMSMixingDepthParameters())


# Initial condition and fluxes
T₀(z) = 20 + Tz*z

models = (cvmix, holtslag, roms, holtslag_roms)

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
            run_until!(model, Δt, (i-1)*10hour)
        end
    end

    @printf("""
        t : %.1f hours

          mixing depths
          =============

               cvmix h : %.2f
            holtslag h : %.2f
                roms h : %.2f
    holtslag + roms  h : %.2f
    \n""", time(cvmix)/hour,
    cvmix.state.h,
    holtslag.state.h,
    roms.state.h,
    holtslag_roms.state.h,
    )

    if i == 1
        vlabel = "CVMix KPP"
        hlabel = "Holtslag \$K\$-profile and CVMix mixing depth model"
        rlabel = "CVMix \$K\$-profile and ROMS mixing depth model"
        mlabel = "Holtslag \$K\$-profile and ROMS mixing depth model"
    else
        vlabel = ""
        hlabel = ""
        rlabel = ""
        mlabel = ""
    end

    if i == 1
        tlabel = text(cvmix.solution.T[end], 0.5,
            @sprintf("\$ t = %.0f \$ hours", time(cvmix)/hour),
            verticalalignment="bottom", horizontalalignment="center", color=defaultcolors[i])
    else
        tlabel = text(maximum(cvmix.solution.T.data)-0.003, -holtslag_roms.state.h,
            @sprintf("\$ t = %.0f \$ hours", time(cvmix)/hour),
            verticalalignment="bottom", horizontalalignment="left", color=defaultcolors[i])
    end

    plot(cvmix.solution.T,          "-",  color=defaultcolors[i], label=vlabel, alpha=0.8, markersize=1.5)
    plot(holtslag.solution.T,       "-.", color=defaultcolors[i], label=hlabel, alpha=0.8, markersize=1.5)
    plot(roms.solution.T,           ":",  color=defaultcolors[i], label=rlabel, alpha=0.8, markersize=1.5)
    plot(holtslag_roms.solution.T,  "--", color=defaultcolors[i], label=mlabel, alpha=0.8, markersize=1.5)
end

legend(fontsize=10)
gcf()
