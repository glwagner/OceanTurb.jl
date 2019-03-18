using
    OceanTurb,
    LinearAlgebra,
    Test

import OceanTurb: Diffusion
import OceanTurb: KPP

#
# Run tests
#

@testset "Utils" begin
    include("utilstests.jl")
    @test test_zeros(Float64)
    @test test_zeros(Float32)
end

@testset "Grids" begin
    include("gridtests.jl")
    nz, Lz = 3, 4.2
    for T in (Float64, Float32, Float16)
        @test test_uniform_grid_type(T, nz, Lz)
        @test test_uniform_grid_spacing(T, nz, Lz)
        @test test_uniform_grid_limits_zc_right(T, nz, Lz)
        @test test_uniform_grid_limits_zf_right(T, nz, Lz)
        @test test_uniform_grid_limits_zc_left(T, nz, Lz)
        @test test_uniform_grid_limits_zf_left(T, nz, Lz)
    end
end

@testset "Fields" begin
    include("fieldtests.jl")
    nz, Lz = 3, 4.2
    for T in (Float64, Float32, Float16)
        @test test_cell_field_construction(T, nz, Lz)
        @test test_face_field_construction(T, nz, Lz)
        @test test_cell_∂z(T)
        @test test_face_∂z(T)
        @test test_cell_plus(T)
        @test test_cell_times(T)
        for loc in (Face, Cell)
            @test test_set_scalar_field(loc, T)
            @test test_set_array_field(loc, T)
            @test test_set_function_field(loc, T)
        end
    end
    @test test_field_indexing()
end

@testset "Diffusion" begin
    include("diffusiontests.jl")
    @test test_diffusion_basic()
    @test test_diffusion_set_c()
    @test test_diffusion_cosine()
end

@testset "KPP" begin
    include("kpptests.jl")
    @test test_constants(g=9.81, α=2.1e-4, ρ₀=1028.1, β=0.99, cP=3904.1, f=1e-4)
    @test test_model_init()
    @test test_Bz()
    @test test_unresolved_KE()
    @test test_surface_layer_average_simple()
    @test test_surface_layer_average_linear()
    @test test_surface_layer_average_steps()
    @test test_Richardson()
    #@test test_mixing_depth()
end
