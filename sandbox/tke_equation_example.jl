using OceanTurb, Printf

function diffusivity!(K, KU, model)
    for i in eachindex(K)
        @inbounds K[i] = KU(model, i)
    end
    return nothing
end

@use_pyplot_utils

usecmbright()

c = Constants(f=1e-4, β=0.0)

modelsetup = (N=128, L=128, stepper=:BackwardEuler, constants=c)

Fb = 0.0 #-1e-7
Fu = -1e-4
Fe = 0.0 #-1e-9
N² = 5e-6
Δt = 1
times = (0, 2, 8, 32) .* hour

K₀ = 1e-5

kpp = ModularKPP.Model(; modelsetup...,
    diffusivity = ModularKPP.LMDDiffusivity(KU₀=K₀, KT₀=K₀, KS₀=K₀)
    )

tke = KPP_TKE.Model(; modelsetup...,
    diffusivity = ModularKPP.LMDDiffusivity(KU₀=K₀, KT₀=K₀, KS₀=K₀),
            tke = KPP_TKE.TKEParameters(Cτ=Inf, CDe=0.5, KU₀=K₀, KT₀=K₀, KS₀=K₀, Ke₀=1e-5)
    )

tkeK = FaceField(tke.grid)
kppK = FaceField(kpp.grid)

# Initial condition and fluxes
Fθ = Fb / (c.α * c.g)
dTdz = N² / (c.α * c.g)
T₀(z) = 20 + dTdz*z
S₀(z) = 20 + dTdz*z

models = (kpp, tke)

for model in models
    model.solution.T = T₀
    model.solution.S = S₀

    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

    model.bcs.U.top = FluxBoundaryCondition(Fu)
end

tke.bcs.e.top = FluxBoundaryCondition(Fe)

fig, axs = subplots(ncols=4, figsize=(8, 6))

sca(axs[1])
removespines("top", "right")
xlabel(L"T")
ylabel(L"z \, \mathrm{(m)}")

sca(axs[2])
removespines("top", "left", "right")
axs[2].tick_params(left=false, labelleft=false)
xlabel(L"U")

#=
sca(axs[3])
removespines("top", "left", "right")
axs[3].tick_params(left=false, labelleft=false)
xlabel(L"V")
=#

sca(axs[3])
removespines("top", "left", "right")
axs[3].tick_params(left=false, labelleft=false)
xlabel(L"e")

sca(axs[4])
removespines("top", "left")
axs[4].tick_params(left=false, labelleft=false, right=true, labelright=true)
xlabel(L"K")
ylabel(L"z \, \mathrm{(m)}")

for i = 1:length(times)
    for model in models
        run_until!(model, Δt, times[i])
    end

    diffusivity!(kppK, ModularKPP.KU, kpp)
    diffusivity!(tkeK, KPP_TKE.K_mixing_time, tke)

    @printf("""
        t : %.1f hours

          surface temperature
          =============

               kpp T(z=0) : %.6f
               tke T(z=0) : %.6f
    \n""", time(kpp)/hour,
    kpp.solution.T[end-1],
    tke.solution.T[end-1],
    )

    if i == 1
        vlabel = (kpp="CVMix KPP", tke="TKE")
    else
        vlabel = (kpp="", tke="")
    end

    #=
    if i == 1
        tlabel = text(kpp.solution.T[end], 0.5,
            @sprintf("\$ t = %.0f \$ hours", time(kpp)/hour),
            verticalalignment="bottom", horizontalalignment="center", color=defaultcolors[i])
    else
        tlabel = text(maximum(kpp.solution.T.data)-0.003, -kpp.state.h,
            @sprintf("\$ t = %.0f \$ hours", time(kpp)/hour),
            verticalalignment="bottom", horizontalalignment="left", color=defaultcolors[i])
    end
    =#

    sca(axs[1])
    plot(kpp.solution.T, "-",  color=defaultcolors[i], label=vlabel.kpp, alpha=0.8, markersize=1.5)
    plot(tke.solution.T, "--", color=defaultcolors[i], label=vlabel.tke, alpha=0.8, markersize=1.5)

    sca(axs[2])
    plot(kpp.solution.U, "-",  color=defaultcolors[i], label=vlabel.kpp, alpha=0.8, markersize=1.5)
    plot(tke.solution.U, "--", color=defaultcolors[i], label=vlabel.tke, alpha=0.8, markersize=1.5)

    #=
    sca(axs[3])
    plot(kpp.solution.V, "-",  color=defaultcolors[i], label=vlabel.kpp, alpha=0.8, markersize=1.5)
    plot(tke.solution.V, "--", color=defaultcolors[i], label=vlabel.tke, alpha=0.8, markersize=1.5)
    =#

    sca(axs[3])
    plot(tke.solution.e, "-", color=defaultcolors[i], label=vlabel.tke, alpha=0.8, markersize=1.5)

    sca(axs[4])
    plot(kppK, "-", color=defaultcolors[i], label=vlabel.kpp, alpha=0.8, markersize=1.5)
    plot(tkeK, "--", color=defaultcolors[i], label=vlabel.tke, alpha=0.8, markersize=1.5)
end

#=
for ax in axs
    sca(ax)
    legend(fontsize=10)
end
=#

gcf()
