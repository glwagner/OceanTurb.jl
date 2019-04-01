#
# Tests for the Diffusion module
#

function test_edmf_basic()
    parameters = EDMF0.Parameters(Cκ=0.42)
    model = EDMF0.Model(N=4, L=2, parameters=parameters)
    model.parameters.Cκ == 0.42
end
