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

    Cεs = [0.01, 0.1, 0.3, 0.5, 1.0]
    avg = zeros(length(Cεs))
    for (j, Cε) in enumerate(Cεs)
        avg[j] = KPP.surface_layer_average(model.solution.T, Cε, i)
    end

    !any(@. !isapprox(avg, T_test))
end

function test_surface_layer_average_linear(;
     N = 10,
     L = 4.1,
     γ = 7.9
    )

    Cε = 1.0 # test only works with full depth layer fraction

    model = KPP.Model(N=N, L=L)

    T_test(z) = γ*z
    model.solution.T = T_test

    avg = FaceField(model.grid)

    for i = 1:N
        avg[i] = KPP.surface_layer_average(model.solution.T, Cε, i)
    end

    h = -model.grid.zf
    avg_answer = @. -0.5 * γ * Cε * h

    @show avg.data avg_answer

    isapprox(avg.data, avg_answer)
end

function test_surface_layer_average_steps()
    Cε = 0.5
     N = 10
     L = 10.0

    model = KPP.Model(N=N, L=L)

    model.solution.T.data[1:8] .= 0
    model.solution.T.data[9] = 2
    model.solution.T.data[10] = 3

    avg_answer = zeros(N+1)
    avg_answer[1:7] .= [10 / (11-i) for i=1:7]
    avg_answer[8] = 4/1.5
    avg_answer[9:10] .= 3
    avg_answer[11] = 0.0

    avg = FaceField(model.grid)

    for i = 1:N
        avg[i] = KPP.surface_layer_average(model.solution.T, Cε, i)
    end

    isapprox(avg.data, avg_answer)
end

function test_Bz(; γ=0.01, g=9.81, ρ₀=1028, α=2e-4, β=0.0, N=10, L=1.0)
    model = KPP.Model(N=N, L=L)
    T₀(z) = γ*z
    model.solution.T = T₀
    Bz_answer = g * α * γ
    Bz = KPP.∂B∂z(g, α, β, model.solution.T, model.solution.S, 3)
    isapprox(Bz, Bz_answer)
end

function test_unresolved_KE(; CKE=0.1, Fb=1e-7, γ=0.01, g=9.81, ρ₀=1028, α=2e-4, β=0.0, N=10, L=1.0)
    Bz = g * α * γ
    i = 2
    T₁(z) = γ*z
    T₂(z) = -γ*z

    model = KPP.Model(N=N, L=L)
    model.solution.T = T₁
    h = -model.grid.zf[i]

    # Test for Bz > 0
    ke₁ = KPP.unresolved_kinetic_energy(CKE, Fb, g, α, β, model.solution.T, model.solution.S, i)
    ke_answer = CKE * h^(4/3) * sqrt(Bz) * Fb^(1/3)

    # Test for Bz < 0
    model.solution.T = T₂
    ke₂ = KPP.unresolved_kinetic_energy(CKE, Fb, g, α, β, model.solution.T, model.solution.S, i)

    ke₂ == 0 && ke₁ ≈ ke_answer
end

function test_Richardson(;
      g = 9.81,
      α = 2.1e-4,
     ρ₀ = 1028.1,
    CRi = 0.69,
    CKE = 1.04,
     Tz = 0.01,
      N = 10,
      L = 20
     )

    parameters = KPP.Parameters(CRi=CRi, CKE=CKE, Cε=0.1/L)
    constants = KPP.Constants(g=g, α=α, ρ₀=ρ₀)

    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)

    T₀(z) = Tz*z
    model.solution.T = T₀

    Ri = FaceField(model.grid)

    @show length(Ri)
    for i = interior(Ri)
        Ri[i] = KPP.Richardson(model, i)
    end

    @show Ri.data

    true
end

function test_mixing_depth(;
      g = 9.81,
      α = 2.1e-4,
     ρ₀ = 1028.1,
    CRi = 0.69,
    CKE = 1.04,
     Tz = 0.01,
      N = 100,
      L = 20
    )

    parameters = KPP.Parameters(CRi=CRi, CKE=CKE, Cε=0.1/L)
    constants = KPP.Constants(g=g, α=α, ρ₀=ρ₀)

    model = KPP.Model(N=N, L=L, parameters=parameters, constants=constants)

    T₀(z) = Tz*z
    model.solution.T = T₀

    h = KPP.mixing_depth(model)
    h_answer = (CRi*CKE)^(3/2) * (α*g*Tz / ρ₀)^(3/4)

    h == h_answer
end
