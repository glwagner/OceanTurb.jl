#
# Tests for the TKE + mass flux module
#

function test_default_model_init(N=4, L=4.3)
    model = TKEMassFlux.Model(grid=UniformGrid(N=N, L=L))
    return model.grid.N == N && model.grid.L == L
end

function time_step_tke_mass_flux_model(boundary_layer_depth, mixing_length, nonlocal_flux, tke_equation, tke_wall_model)

    model = TKEMassFlux.Model(                grid = UniformGrid(N=4, L=3), 
                              boundary_layer_depth = boundary_layer_depth,
                                     mixing_length = mixing_length,
                                     nonlocal_flux = nonlocal_flux,
                                      tke_equation = tke_equation,
                                    tke_wall_model = tke_wall_model,
                                           stepper = :BackwardEuler)
    iterate!(model, 1e-16) 
    return true
end

boundary_layer_depth_models = (
    nothing,
    ModularKPP.LMDMixingDepth(),
    ModularKPP.ROMSMixingDepth()
)

nonlocal_flux_models = (
    nothing,
    ModularKPP.LMDCounterGradientFlux(),
    ModularKPP.DiagnosticPlumeModel()
)

tke_equation_models = (
    TKEMassFlux.TKEParameters(),
)

mixing_length_models = (
    TKEMassFlux.SimpleMixingLength(),
    TKEMassFlux.TanEtAl2018MixingLength(),
    TKEMassFlux.EquilibriumMixingLength(),
)

wall_models = (
    nothing,
    TKEMassFlux.PrescribedNearWallTKE(),
    TKEMassFlux.SurfaceValueScaling(),
    TKEMassFlux.SurfaceFluxScaling(),
)

@testset "TKEMassFlux" begin
    @test test_default_model_init()

    # Modularity
    for boundary_layer_depth in boundary_layer_depth_models
        for nonlocal_flux in nonlocal_flux_models
            for tke_equation in tke_equation_models
                for mixing_length in mixing_length_models
                    for wall_model in wall_models
                        @test time_step_tke_mass_flux_model(boundary_layer_depth, mixing_length, 
                                                            nonlocal_flux, tke_equation, wall_model)
                    end
                end
            end
        end
    end
end
