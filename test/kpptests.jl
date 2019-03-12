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

function test_mixing_depth(;
      g = 9.81,
      α = 2.1e-4,
     ρ₀ = 1028.1,
    CRi = 0.69,
    CKE = 1.04,
     Tz = 0.01
    )

    parameters = KPP.Parameters(CRi=CRi, CKE=CKE, Cε=0.0001)
    constants = KPP.Constants(g=g, α=α, ρ₀=ρ₀)

    model = KPP.Model(N=100, L=100.0, parameters=parameters)

    T(z) = Tz*z
    model.solution.T = T

    h = KPP.mixing_depth(model)

    @show h h/model.grid.L
    
    h == (CRi*CKE)^(3/2) * (α*g*Tz / ρ₀)^(3/4)
end


