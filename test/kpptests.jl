#
# Tests for the K-Profile-Parameterization
#

function test_constants(; kwargs...)
    constants = KPP.Constants(; kwargs...)
    result = true
    for (k, v) in kwargs
        if !result
            result = getproperty(constants, k) == v
        end
    end
    return result
end

function test_model_init(N=4, L=4.3)
    model = KPP.Model(N=N, L=L)
    model.grid.N == N && model.grid.L == L
end

function test_surface_layer_average_simple(;
    N = 10,
    L = 1.1,
    i = 9
    )

    model = KPP.Model(N=N, L=L)

    T_test = 2.1
    model.solution.T = T_test

    CSLs = [0.01, 0.1, 0.3, 0.5, 1.0]
    avg = zeros(length(CSLs))
    for (j, CSL) in enumerate(CSLs)
        avg[j] = KPP.surface_layer_average(model.solution.T, CSL, i)
    end

    !any(@. !isapprox(avg, T_test))
end

function test_surface_layer_average_linear(; N=10, L=4.1, γ=7.9)
    CSL = 1.0 # test only works with full depth layer fraction
    model = KPP.Model(N=N, L=L)
    T_test(z) = γ*z
    model.solution.T = T_test

    avg = FaceField(model.grid)

    for i = 1:N
        avg[i] = KPP.surface_layer_average(model.solution.T, CSL, i)
    end

    h = -model.grid.zf
    avg_answer = @. -0.5 * γ * CSL * h
    isapprox(avg.data, avg_answer)
end

function analytical_surface_layer_average()
    CSL, N, L = 0.5, 10, 10.0
    T, avg = zeros(N), zeros(N+1)

    T[1:8] .= 0
    T[9] = 2
    T[10] = 3

    avg[1:7] .= [10 / (11-i) for i=1:7]
    avg[8] = 4/1.5
    avg[9:10] .= 3
    avg[11] = 0.0

    return CSL, N, L, T, avg
end

function test_surface_layer_average_steps()
    CSL, N, L, T, avg_answer = analytical_surface_layer_average()
    model = KPP.Model(N=N, L=L)
    model.solution.T = T

    avg = zeros(N+1)
    for i = 1:N
        avg[i] = KPP.surface_layer_average(model.solution.T, CSL, i)
    end

    avg ≈ avg_answer
end

function test_surface_layer_average_epsilon(; γ=2.1, CSL=0.01, N=10, L=1.0, T₀₀=281.2)
    T₀(z) = T₀₀*γ*z
    model = KPP.Model(N=N, L=L)
    model.solution.T = T₀

    i = 2
    h = -model.grid.zf[i]
    avg = KPP.surface_layer_average(model.solution.T, CSL, i)
    avg_answer = T₀(-model.grid.Δf/2)

    avg ≈ avg_answer
end

function test_Δ0()
    model = KPP.Model()
    (KPP.Δ(model.solution.U, 0.01, 2) == 0 &&
     KPP.Δ(model.solution.V, 0.01, 2) == 0 &&
     KPP.Δ(model.solution.T, 0.01, 2) == 0 &&
     KPP.Δ(model.solution.S, 0.01, 2) == 0)
end


function test_Δ1()
    CSL, N, L, T, surf_avg = analytical_surface_layer_average()
    parameters = KPP.Parameters(CSL=CSL)
    model = KPP.Model(N=N, L=L, parameters=parameters)
    model.solution.T = T

    Δ_lower = surf_avg[2]
    Δ_upper = surf_avg[10] - 2.5

    KPP.Δ(model.solution.T, CSL, 2) == Δ_lower && KPP.Δ(model.solution.T, CSL, 10) == Δ_upper
end

function test_Δ2(; γ=2.1, CSL=0.001, N=10, L=1.0, i=N)
    model = KPP.Model(N=N, L=L)
    T₀(z) = γ*z
    model.solution.T = T₀
    T = model.solution.T
    h = -model.grid.zf[i]
    Δ_answer = T[N] - T₀(-h)
    KPP.Δ(model.solution.T, CSL, i) ≈ Δ_answer
end


function test_Δ3(; CSL=0.5, N=20, L=20)
    U₀ = 3
    parameters = KPP.Parameters(CSL=CSL)
    model = KPP.Model(N=N, L=L, parameters=parameters)
    U, V, T, S = model.solution

    ih = 12
    @views U.data[ih:N] .= U₀
    @views U.data[1:ih-1] .= -U₀

    KPP.Δ(U, CSL, ih) == U₀
end


function test_buoyancy_gradient(; γ=0.01, g=9.81, ρ₀=1028, α=2e-4, β=0.0, N=10, L=1.0)
    model = KPP.Model(N=N, L=L)
    T₀(z) = γ*z
    model.solution.T = T₀
    Bz_answer = g * α * γ
    Bz = KPP.∂B∂z(model.solution.T, model.solution.S, g, α, β, 3)
    Bz ≈ Bz_answer
end

function test_unresolved_KE(; CKE=0.1, CKE₀=1e-11, Fb=1e-7, γ=0.01, g=9.81,
                            ρ₀=1028, α=2e-4, β=0.0, N=10, L=1.0)
    Bz = g * α * γ
    T₁(z) = γ*z
    T₂(z) = -γ*z
    model = KPP.Model(N=N, L=L)

    i = 2
    h = -model.grid.zf[i]

    # Test for Bz > 0
    model.solution.T = T₁
    Bz = KPP.∂B∂z(model, i)
    ke₁ = KPP.unresolved_kinetic_energy(h, Bz, Fb, CKE, CKE₀, g, α, β, i)
    ke₁_answer = CKE * h^(4/3) * sqrt(Bz) * Fb^(1/3) + CKE₀

    # Test for Bz < 0
    model.solution.T = T₂
    Bz = KPP.∂B∂z(model, i)
    ke₂ = KPP.unresolved_kinetic_energy(h, Bz, Fb, CKE, CKE₀, g, α, β, i)
    ke₂_answer = CKE₀

    ke₂ ≈ ke₂_answer && ke₁ ≈ ke₁_answer
end

function test_update_state(; N=10, L=20, Fθ=5.1e-3)
    model = KPP.Model(N=N, L=L)
    temperature_bc = FluxBoundaryCondition(Fθ)
    model.bcs.T.top = temperature_bc

    KPP.update_state!(model)

    model.state.Fθ == Fθ
end

function test_bulk_richardson_number(; g=9.81, α=2.1e-4, CRi=0.3, CKE=1.04,
                                       CKE₀=0.0, γ=0.01, N=20, L=20, Fb=2.1e-5)
    parameters = KPP.Parameters(CRi=CRi, CKE=CKE, CKE₀=CKE₀, CSL=0.1/N)
    constants = KPP.Constants(g=g, α=α)
    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)
    U, V, T, S = model.solution

    # Initial condition and top flux
    T₀(z) = γ*z
    model.solution.T = T₀

    # "Surface" temperature and buoyancy
    T_N = T[N]
    B_N = g*α*T[N]

    Fθ = Fb / (g*α)
    top_bc_T = FluxBoundaryCondition(Fθ)
    model.bcs.T.top = top_bc_T
    KPP.update_state!(model)

    Ri = FaceField(model.grid)
    Ri_answer = FaceField(model.grid)
    for i = interiorindices(Ri)
        Ri[i] = KPP.bulk_richardson_number(model, i)

        h = - model.grid.zf[i]
        h⁺ = h * (1 - 0.5*model.parameters.CSL)
        Bz = KPP.∂B∂z(model, i)
        uke = KPP.unresolved_kinetic_energy(h, Bz, Fb, CKE, CKE₀, g, α, model.constants.β, i)
        Ri_answer[i] = h⁺ * (B_N - g*α*T₀(-h)) / uke
    end

    isapprox(Ri_answer.data, Ri.data)
end

function test_mixing_depth_convection(; g=9.81, α=2.1e-4, CRi=0.3, CKE=1.04,
                                       γ=0.01, N=200, L=20, Fb=4.1e-5)
    parameters = KPP.Parameters(CRi=CRi, CKE=CKE, CKE₀=0.0, CSL=0.1/N)
    constants = KPP.Constants(g=g, α=α)
    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)
    U, V, T, S = model.solution

    T₀(z) = γ*z
    model.solution.T = T₀
    model.solution.T[N] = 0

    Fθ = Fb / (g*α)
    temperature_bc = FluxBoundaryCondition(Fθ)
    model.bcs.T.top = temperature_bc
    KPP.update_state!(model)

    Ri = FaceField(model.grid)
    for i = interiorindices(Ri)
        Ri[i] = KPP.bulk_richardson_number(model, i)
    end

    h = KPP.mixing_depth(model)
    h_answer = (CRi*CKE)^(3/2) * (α*g*γ)^(-3/4) * sqrt(Fb) / (1 - 0.5*model.parameters.CSL)

    isapprox(h, h_answer, rtol=1e-3)
end

function test_mixing_depth_shear(; CSL=0.5, N=20, L=20, CRi=1.0)
    T₀ = 1
    U₀ = 3
    parameters = KPP.Parameters(CRi=CRi, CSL=CSL, CKE₀=0.0)
    constants = KPP.Constants(g=1, α=1)
    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)
    U, V, T, S = model.solution

    h = CRi * U₀^2 / T₀ / (1 - 0.5*model.parameters.CSL)
    ih = Int(N*(1 - h/L) + 1)
    @views T.data[ih:N] .= T₀
    @views U.data[ih:N] .= U₀

    @views T.data[1:ih-1] .= -T₀
    @views U.data[1:ih-1] .= -U₀

    KPP.mixing_depth(model) ≈ h
end

function test_unstable()
    model = KPP.Model()

    model.bcs.T.top = FluxBoundaryCondition(1)
    KPP.update_state!(model)
    stability₊ = KPP.isunstable(model)

    model.bcs.T.top = FluxBoundaryCondition(-1)
    KPP.update_state!(model)
    stability₋ = KPP.isunstable(model)

    stability₊ == true && stability₋ == false
end

function test_zero_turbulent_velocity()
    model = KPP.Model()
    KPP.ωτ(model) == 0.0 &&  KPP.ωb(model) == 0.0
end

function test_friction_velocity()
    model = KPP.Model()
    model.bcs.U.top = FluxBoundaryCondition(sqrt(8))
    model.bcs.V.top = FluxBoundaryCondition(-sqrt(8))
    KPP.update_state!(model)
    KPP.ωτ(model) ≈ 2
end

function test_convective_velocity()
    model = KPP.Model()
    γ = 0.01
    T₀(z) = γ*z
    model.solution.T = T₀

    Fb = 2.1
    Fθ = Fb / (model.constants.α*model.constants.g)
    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    KPP.update_state!(model)

    h = KPP.mixing_depth(model)

    KPP.ωb(model) ≈ (h*Fb)^(1/3)
end

function test_turb_velocity_pure_convection(N=20, L=20, Cb_U=3.1, Cb_T=1.7, CSL=1e-16)
    # Zero wind + convection => w_scale_U = Cb_U * CSL^(1/3) * ωb.
    parameters = KPP.Parameters(CRi=1.0, CKE=1.0, CKE₀=0.0, CSL=CSL, Cb_U=Cb_U, Cb_T=Cb_T)
    constants = KPP.Constants(g=1, α=1)
    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)
    U, V, T, S = model.solution

    T₀(z) = z
    model.solution.T = T₀
    model.solution.T[N] = 0

    Fb = 100
    model.bcs.T.top = FluxBoundaryCondition(Fb)
    KPP.update_state!(model)

    i = 16
    h = sqrt(Fb) / (1-0.5CSL) # requires h to be an integer... ?
    ωb = (h*Fb)^(1/3)

    (KPP.w_scale_U(model, i) ≈ Cb_U * CSL^(1/3) * ωb &&
     KPP.w_scale_V(model, i) ≈ Cb_U * CSL^(1/3) * ωb &&
     KPP.w_scale_T(model, i) ≈ Cb_T * CSL^(1/3) * ωb &&
     KPP.w_scale_S(model, i) ≈ Cb_T * CSL^(1/3) * ωb )
end

function test_turb_velocity_pure_wind(; CSL=0.5, Cτ=0.7, N=20, L=20, CRi=1.0)
    T₀ = 1
    U₀ = 3
    parameters = KPP.Parameters(CRi=CRi, Cτ=Cτ, CSL=CSL)
    constants = KPP.Constants(g=1, α=1)
    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)
    U, V, T, S = model.solution

    h = CRi * U₀^2 / T₀
    ih = Int(N*(1 - h/L) + 1)
    @views T.data[ih:N] .= T₀
    @views U.data[ih:N] .= U₀

    @views T.data[1:ih-1] .= -T₀
    @views U.data[1:ih-1] .= -U₀

    Fu = 2.1
    model.bcs.U.top = FluxBoundaryCondition(Fu)
    KPP.update_state!(model)

    w_scale = Cτ * sqrt(Fu)

    (KPP.w_scale_U(model, 3) == w_scale &&
     KPP.w_scale_V(model, 3) == w_scale &&
     KPP.w_scale_T(model, 3) == w_scale &&
     KPP.w_scale_S(model, 3) == w_scale )
end


function test_turb_velocity_wind_stab(; CSL=0.5, Cτ=0.7, N=20, L=20, CRi=1.0, Cstab=0.3)
    T₀ = 1
    U₀ = 3
    parameters = KPP.Parameters(CRi=CRi, Cτ=Cτ, CSL=CSL, Cstab=Cstab)
    constants = KPP.Constants(g=1, α=1)
    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)
    U, V, T, S = model.solution

    h = CRi * U₀^2 / T₀
    ih = Int(N*(1 - h/L) + 1) # 12
    @views T.data[ih:N] .= T₀
    @views U.data[ih:N] .= U₀

    @views T.data[1:ih-1] .= -T₀
    @views U.data[1:ih-1] .= -U₀

    Fu = 2.1
    model.bcs.U.top = FluxBoundaryCondition(Fu)
    Fθ = -1.3
    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    KPP.update_state!(model)

    rb = abs(h*Fθ) / Fu^(3/2)

    id = 16 # d=5/9
    d = 5/9
    w_scale = Cτ * sqrt(Fu) / (1 + Cstab * rb * d)

    (KPP.w_scale_U(model, id) ≈ w_scale &&
     KPP.w_scale_V(model, id) ≈ w_scale &&
     KPP.w_scale_T(model, id) ≈ w_scale &&
     KPP.w_scale_S(model, id) ≈ w_scale )
end

function test_turb_velocity_wind_unstab(; CKE=0.0, CSL=0.5, Cτ=0.7, N=20,
                                        L=20, CRi=(1-0.5CSL), Cunst=0.3)
    T₀ = 1
    U₀ = 3
    parameters = KPP.Parameters(CRi=CRi, Cτ=Cτ, CKE=CKE, CKE₀=0.0,
                                CSL=CSL, Cunst=Cunst, Cd_U=Inf, Cd_T=Inf)
    constants = KPP.Constants(g=1, α=1)
    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)
    U, V, T, S = model.solution

    h = CRi * U₀^2 / T₀ / (1-0.5CSL)
    ih = Int(N*(1 - h/L) + 1) # 12
    @views T.data[ih:N] .= T₀
    @views U.data[ih:N] .= U₀

    @views T.data[1:ih-1] .= -T₀
    @views U.data[1:ih-1] .= -U₀

    Fu = 2.1
    model.bcs.U.top = FluxBoundaryCondition(Fu)
    Fθ = 1.1
    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    KPP.update_state!(model)

    rb = abs(h*Fθ) / (Fu)^(3/2)
    id1 = 16 # d=5/9
    id2 = 18 # d=3/9
    d1 = 5/9
    d2 = 3/9

    w_scale_U1 = Cτ * sqrt(Fu) * (1 + Cunst * rb * CSL)^(1/4)
    w_scale_T1 = Cτ * sqrt(Fu) * (1 + Cunst * rb * CSL)^(1/2)

    w_scale_U2 = Cτ * sqrt(Fu) * (1 + Cunst * rb * d2)^(1/4)
    w_scale_T2 = Cτ * sqrt(Fu) * (1 + Cunst * rb * d2)^(1/2)

    (KPP.w_scale_U(model, id1) ≈ w_scale_U1 &&
     KPP.w_scale_V(model, id1) ≈ w_scale_U1 &&
     KPP.w_scale_T(model, id1) ≈ w_scale_T1 &&
     KPP.w_scale_S(model, id1) ≈ w_scale_T1 &&
     KPP.w_scale_U(model, id2) ≈ w_scale_U2 &&
     KPP.w_scale_V(model, id2) ≈ w_scale_U2 &&
     KPP.w_scale_T(model, id2) ≈ w_scale_T2 &&
     KPP.w_scale_S(model, id2) ≈ w_scale_T2 )
end

function test_conv_velocity_wind(; CKE=0.0, CKE₀=0.0, CSL=0.5, Cτ=0.7, N=20, L=20, CRi=(1-0.5CSL),
                                 Cb_U=1.1, Cb_T=0.1)
    T₀ = 1
    U₀ = 3

    parameters = KPP.Parameters(CRi=CRi, Cτ=Cτ, CKE=CKE, CKE₀=0.0, CSL=CSL, Cd_U=0.0, Cd_T=0.0,
                                Cb_U=1.1, Cb_T=0.1)

    constants = KPP.Constants(g=1, α=1)
    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)
    U, V, T, S = model.solution

    h = CRi * U₀^2 / T₀ / (1-0.5CSL)
    ih = Int(N*(1 - h/L) + 1) # 12
    @views T.data[ih:N] .= T₀
    @views U.data[ih:N] .= U₀

    @views T.data[1:ih-1] .= -T₀
    @views U.data[1:ih-1] .= -U₀

    Fu = -2.1
    model.bcs.U.top = FluxBoundaryCondition(Fu)
    Fθ = 0.5
    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    KPP.update_state!(model)

    rb = abs(h*Fθ) / abs(Fu)^(3/2)
    rτ = 1/rb
    id1 = 16 # d=5/9
    id2 = 18 # d=3/9
    d1 = 5/9
    d2 = 3/9

    Cτb_U = model.parameters.Cτb_U
    Cτb_T = model.parameters.Cτb_T

    w_scale_U1 = Cb_U * abs(h*Fθ)^(1/3) * (CSL + Cτb_U * rτ)^(1/3)
    w_scale_T1 = Cb_T * abs(h*Fθ)^(1/3) * (CSL + Cτb_T * rτ)^(1/3)

    w_scale_U2 = Cb_U * abs(h*Fθ)^(1/3) * (d2 + Cτb_U * rτ)^(1/3)
    w_scale_T2 = Cb_T * abs(h*Fθ)^(1/3) * (d2 + Cτb_T * rτ)^(1/3)

    (KPP.w_scale_U(model, id1) ≈ w_scale_U1 &&
     KPP.w_scale_V(model, id1) ≈ w_scale_U1 &&
     KPP.w_scale_T(model, id1) ≈ w_scale_T1 &&
     KPP.w_scale_S(model, id1) ≈ w_scale_T1 &&
     KPP.w_scale_U(model, id2) ≈ w_scale_U2 &&
     KPP.w_scale_V(model, id2) ≈ w_scale_U2 &&
     KPP.w_scale_T(model, id2) ≈ w_scale_T2 &&
     KPP.w_scale_S(model, id2) ≈ w_scale_T2 )
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
    bcs.bottom = FluxBoundaryCondition(bottom_flux)

    c = getproperty(model.solution, fieldname)
    C₀ = integral(c)
    C(t) = C₀ - (top_flux - bottom_flux) * t

    dt = 1e-6
    iterate!(model, dt, 10)

    return C(time(model)) ≈ integral(c)
end
