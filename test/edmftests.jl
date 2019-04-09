#
# Tests for the Diffusion module
#

function test_edmf_basic()
    parameters = EDMF.Parameters(Cκ=0.42)
    model = EDMF.ZeroPlumeModel(N=4, L=2, parameters=parameters)
    model.parameters.Cκ == 0.42
end
