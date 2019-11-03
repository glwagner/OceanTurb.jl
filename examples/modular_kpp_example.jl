using OceanTurb, Printf

using OceanTurb.ModularKPP: HoltslagDiffusivity, ROMSMixingDepth, LMDCounterGradientFlux,
                            DiagnosticPlumeModel, LMDDiffusivity, LMDMixingDepth

@use_pyplot_utils

usecmbright()

modelsetup = (N=256, L=256, stepper=:BackwardEuler)

Fb = 1e-7
Fu = 0.0
N² = 1e-5
Δt = 2minute
t_plot = (0hour, 12hour, 48hour)

name = "Free convection
    \\small{with \$ \\overline{w b} |_{z=0} = 10^{-7} \\, \\mathrm{m^2 \\, s^{-3}}\$}"

# The standard setup.
        cvmix = ModularKPP.Model(; modelsetup...,
                                    diffusivity = LMDDiffusivity(),
                                    mixingdepth = LMDMixingDepth(),
                                   nonlocalflux = LMDCounterGradientFlux())

# The standard setup except with a plume model rather than a counter-gradient flux model.
 cvmix_plumes = ModularKPP.Model(; modelsetup...,
                                    diffusivity = LMDDiffusivity(),
                                    mixingdepth = LMDMixingDepth(),
                                   nonlocalflux = DiagnosticPlumeModel())

# Standard nonlocal flux and mixing depth model, but a much simpler diffusivity model.
     holtslag = ModularKPP.Model(; modelsetup...,
                                    diffusivity = HoltslagDiffusivity(),
                                    mixingdepth = LMDMixingDepth(),
                                   nonlocalflux = LMDCounterGradientFlux())

# Use the mixing depth model employed by ROMS, the Regional Ocean Modeling System.
         roms = ModularKPP.Model(; modelsetup...,
                                    diffusivity = LMDDiffusivity(),
                                    mixingdepth = ROMSMixingDepth(),
                                   nonlocalflux = LMDCounterGradientFlux())

# Pair the simple Holtslag diffusivity model with the ROMS mixing depth model.
holtslag_roms = ModularKPP.Model(; modelsetup...,
                                    diffusivity = HoltslagDiffusivity(),
                                    mixingdepth = ROMSMixingDepth(),
                                   nonlocalflux = LMDCounterGradientFlux())


# Initial condition and fluxes
dTdz = N² / (cvmix.constants.α * cvmix.constants.g)
T₀(z) = 20 + dTdz*z

models = (cvmix, holtslag, roms, holtslag_roms, cvmix_plumes)

for model in models
    model.solution.T = T₀

    Fθ = Fb / (model.constants.α * model.constants.g)
    model.bcs.U.top = FluxBoundaryCondition(Fu)

    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

end

fig, axs = subplots()

removespines("top", "right")
xlabel("Temperature \$ \\, {}^\\circ \\mathrm{C} \$")
ylabel(L"z \, \mathrm{(m)}")

for i = 1:length(t_plot)
    for model in models
        run_until!(model, Δt, t_plot[i])
    end

    @printf("""
        t : %.1f hours

          mixing depths
          =============

               cvmix h : %.2f
      cvmix + plumes h : %.2f
            holtslag h : %.2f
                roms h : %.2f
    holtslag + roms  h : %.2f
    \n""", time(cvmix)/hour,
    cvmix.state.h,
    cvmix_plumes.state.h,
    holtslag.state.h,
    roms.state.h,
    holtslag_roms.state.h,
    )

    if i == 1
        labels = [
            "CVMix KPP",
            "CVMix KPP with a plume model",
            "Holtslag \$K\$-profile and CVMix mixing depth model",
            "CVMix \$K\$-profile and ROMS mixing depth model",
            "Holtslag \$K\$-profile and ROMS mixing depth model"
           ]
    else
        labels = ["" for i in 1:5]
    end

    if i == 1
        tlabel = text(cvmix.solution.T[end], 0.5,
            @sprintf("\$ t = %.0f \$ hours", time(cvmix)/hour),
            verticalalignment="bottom", horizontalalignment="center", color=defaultcolors[i])
    else
        tlabel = text(maximum(cvmix.solution.T.data)-0.001, -holtslag_roms.state.h,
            @sprintf("\$ t = %.0f \$ hours", time(cvmix)/hour),
            verticalalignment="bottom", horizontalalignment="left", color=defaultcolors[i])
    end

    plot(cvmix.solution.T,           "-", color=defaultcolors[i], label=labels[1], alpha=0.8, markersize=1.5)
    plot(cvmix_plumes.solution.T,   "--", color=defaultcolors[i], label=labels[2], alpha=0.8, markersize=1.5)
    plot(holtslag.solution.T,       "-.", color=defaultcolors[i], label=labels[3], alpha=0.8, markersize=1.5)
    plot(roms.solution.T,            "^", color=defaultcolors[i], label=labels[4], alpha=0.8, markersize=1.5)
    plot(holtslag_roms.solution.T,   ":", color=defaultcolors[i], label=labels[5], alpha=0.8, markersize=1.5)
end

title(name)
legend(fontsize=10)
gcf()

#savefig("figs/free_convection_intermodel.png", dpi=480)
