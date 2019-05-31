#
# Tests for the moular K-Profile-Parameterization module
#

function test_model_init(N=4, L=4.3)
    model = ModularKPP.Model(N=N, L=L)
    model.grid.N == N && model.grid.L == L
end

function test_diffusivity_plain(; K₀=1.1)
    parameters = KPP.Parameters(KU₀=K₀, KT₀=K₀, KS₀=K₀)
    model = KPP.Model(parameters=parameters)
    KPP.update_state!(model)

    m = model
    KU = FaceField(model.grid)
    KV = FaceField(model.grid)
    KT = FaceField(model.grid)
    KS = FaceField(model.grid)
    for i = 1:model.grid.N+1
        KU[i] = KPP.KU(model, i)
        KV[i] = KPP.KV(model, i)
        KT[i] = KPP.KT(model, i)
        KS[i] = KPP.KS(model, i)
    end

    (!any(@. KU.data != K₀) &&
     !any(@. KV.data != K₀) &&
     !any(@. KT.data != K₀) &&
     !any(@. KS.data != K₀) )
end


function test_kpp_diffusion_cosine(stepper=:ForwardEuler)
    parameters = KPP.Parameters(KT₀=1.0, KS₀=1.0)
    model = KPP.Model(N=100, L=π/2, parameters=parameters, stepper=stepper)
    z = model.grid.zc

    c_init(z) = cos(2z)
    c_ans(z, t) = exp(-4t) * c_init(z)

    model.solution.T = c_init
    model.solution.S = c_init

    KPP.update_state!(model)
    U, V, T, S = model.solution

    m = model
    i = 3
    dt = 1e-3
    iterate!(model, dt)

    # The error tolerance is a bit arbitrary.
    norm(c_ans.(z, time(model)) .- data(model.solution.T)) < model.grid.N*1e-6
end

function test_flux(stepper=:ForwardEuler; fieldname=:U, top_flux=0.3, bottom_flux=0.13, N=10)
    model = KPP.Model(N=N, L=1, stepper=stepper)

    bcs = getproperty(model.bcs, fieldname)
    bcs.top = FluxBoundaryCondition(top_flux)
    bottom_K = getproperty(model.timestepper.eqn.K, fieldname)(model, 1)
    bcs.bottom = GradientBoundaryCondition(-bottom_flux / bottom_K)

    c = getproperty(model.solution, fieldname)
    C₀ = integral(c)
    C(t) = C₀ - (top_flux - bottom_flux) * t

    dt = 1e-6
    iterate!(model, dt, 10)

    return C(time(model)) ≈ integral(c)
end
