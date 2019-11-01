#
# Tests for the moular K-Profile-Parameterization module
#

function test_model_init(N=4, L=4.3)
    model = ModularKPP.Model(N=N, L=L)
    model.grid.N == N && model.grid.L == L
end

@testset "Modular KPP" begin
    @test test_model_init()
end
