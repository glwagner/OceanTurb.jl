#
# Tests for the TKE + mass flux module
#

function test_default_model_init(N=4, H=4.3)
    model = TKEMassFlux.Model(grid=UniformGrid(N=N, H=H))
    return model.grid.N == N && model.grid.H == H
end

function time_step_tke_mass_flux_model(boundary_layer_depth, mixing_length, nonlocal_flux,
                                       tke_equation, tke_wall_model, eddy_diffusivities)

    model = TKEMassFlux.Model(                grid = UniformGrid(N=4, H=4), 
                              boundary_layer_depth = boundary_layer_depth,
                                     mixing_length = mixing_length,
                                     nonlocal_flux = nonlocal_flux,
                                eddy_diffusivities = eddy_diffusivities,
                                      tke_equation = tke_equation,
                                    tke_wall_model = tke_wall_model,
                                           stepper = :BackwardEuler)
    time_step!(model, 1e-16, 3)

    return true
end

eddy_diffusivity_models = (
                           TKEMassFlux.SinglePrandtlDiffusivities(),
                           TKEMassFlux.IndependentDiffusivities(),
                           TKEMassFlux.RiDependentDiffusivities(),
                          )

boundary_layer_depth_models = (
                               nothing,
                               ModularKPP.LMDMixingDepth(),
                               ModularKPP.ROMSMixingDepth(),
                              )

nonlocal_flux_models = (
                        nothing,
                        ModularKPP.LMDCounterGradientFlux(),
                        ModularKPP.DiagnosticPlumeModel(),
                       )

tke_equation_models = (
                       TKEMassFlux.TKEParameters(),
                      )

mixing_length_models = (
                        TKEMassFlux.SimpleMixingLength(),
                        TKEMassFlux.EquilibriumMixingLength(),
                       )

wall_models = (
               nothing,
               TKEMassFlux.PrescribedNearWallTKE(),
               TKEMassFlux.PrescribedSurfaceTKEValue(),
               TKEMassFlux.PrescribedSurfaceTKEFlux(),
              )

@testset "TKEMassFlux" begin
    @test test_default_model_init()

    println("TKEMassFlux tests begin:")
    println(" ")
    
    # Modularity
    for boundary_layer_depth in boundary_layer_depth_models
        for nonlocal_flux in nonlocal_flux_models
            for tke_equation in tke_equation_models
                for mixing_length in mixing_length_models
                    for wall_model in wall_models
                        for eddy_diffusivity_model in eddy_diffusivity_models

                            @test time_step_tke_mass_flux_model(boundary_layer_depth,
                                                                mixing_length, 
                                                                nonlocal_flux,
                                                                tke_equation,
                                                                wall_model,
                                                                eddy_diffusivity_model)
                            println("testing... OK") # this is a hack to avoid travis-ci timeout
                        end
                    end
                end
            end
        end
    end
end
