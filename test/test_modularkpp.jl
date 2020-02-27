#
# Tests for the moular K-Profile-Parameterization module
#

function test_default_model_init(N=4, H=4.3)
    model = ModularKPP.Model(grid=UniformGrid(N=N, H=H))
    return model.grid.N == N && model.grid.H == H
end

function time_step_model(diffusivity, nonlocalflux, mixingdepth, kprofile)

    model = ModularKPP.Model(        grid = UniformGrid(N=4, H=4), 
                              diffusivity = diffusivity,
                             nonlocalflux = nonlocalflux,
                              mixingdepth = mixingdepth,
                                 kprofile = kprofile)

    time_step!(model, 1e-16)
    return true
end

diffusivity_models = (
    ModularKPP.LMDDiffusivity(),
    ModularKPP.HoltslagDiffusivity()
)

nonlocal_flux_models = (
    ModularKPP.LMDCounterGradientFlux(),
    ModularKPP.DiagnosticPlumeModel()
)

mixing_depth_models = (
    ModularKPP.LMDMixingDepth(),
    ModularKPP.ROMSMixingDepth()
)

shape_functions = (
    ModularKPP.StandardCubicPolynomial(),
    ModularKPP.GeneralizedCubicPolynomial()
)

@testset "Modular KPP" begin
    @test test_default_model_init()

    # Modularity
    for diffusivity in diffusivity_models
        for nonlocalflux in nonlocal_flux_models
            for mixingdepth in mixing_depth_models
                for kprofile in shape_functions
                    @test time_step_model(diffusivity, nonlocalflux, mixingdepth, kprofile)
                end
            end
        end
    end
end
