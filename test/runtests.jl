using
  OceanTurb,
  Test

import OceanTurb
import OceanTurb: Diffusion
  
# --
# Define tests
# --

function test_zeros(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  OceanTurb.@zeros T dims a2 b2
  eltype(a1) == T && a1 == a2 && b1 == b2
end

function test_uniform_grid_type(T, nz, Lz)
  grid = UniformGrid(T, nz, Lz)
  T == typeof(grid.Lz) && T == eltype(grid.zc) && T == eltype(grid.zf)
end

function test_uniform_grid_spacing(T, nz, Lz)
  dz = T(Lz/nz)
  grid = UniformGrid(T, nz, Lz)
  dz == grid.dzc && dz == grid.dzf
end

function test_uniform_grid_limits_zc_left(T, nz, Lz)
  dz = Lz/nz
  grid = UniformGrid(T, nz, Lz)
  grid.zc[1] ≈ -Lz+0.5dz
end

function test_uniform_grid_limits_zc_right(T, nz, Lz)
  dz = Lz/nz
  grid = UniformGrid(T, nz, Lz)
  grid.zc[nz] ≈ -0.5dz
end

function test_uniform_grid_limits_zf_left(T, nz, Lz)
  dz = Lz/nz
  grid = UniformGrid(T, nz, Lz)
  grid.zf[1] ≈ -Lz
end

function test_uniform_grid_limits_zf_right(T, nz, Lz)
  dz = Lz/nz
  grid = UniformGrid(T, nz, Lz)
  grid.zf[nz+1] ≈ 0
end

function test_cell_field_construction(T, nz, Lz)
  grid = UniformGrid(T, nz, Lz)
  c = CellField(grid)
  size(c.data) == (OceanTurb.cell_field_size(nz),)
end

function test_face_field_construction(T, nz, Lz)
  grid = UniformGrid(T, nz, Lz)
  f = FaceField(grid)
  size(f.data) == (OceanTurb.face_field_size(nz),)
end

function test_field_indexing()
  nz = 3
  Lz = 4.2
  val = 2.1
  grid = UniformGrid(nz, Lz)
  c = CellField(grid)
  f = FaceField(grid)
  c[2] = val
  f[2] = val
  c[2] == val && f[2] == val
end

function test_cell_∂z(T)
  nz = 2
  Lz = 2.0
  grid = UniformGrid(T, nz, Lz)
  c = CellField([2, 4, 6, 8], grid)
  f = ∂z(c)
  f.data == Vector{Float64}([2, 2, 2])
end

function test_face_∂z(T)
  nz = 2
  Lz = 2.0
  grid = UniformGrid(T, nz, Lz)
  f = FaceField([2, 4, 6], grid)
  c = ∂z(f)
  c_answer = CellField([0, 2, 2, 0], grid)
  c.data == c_answer.data
end

function test_cell_plus(T)
  nz = 2
  Lz = 2.0
  grid = UniformGrid(T, nz, Lz)
  c1 = CellField([0, 1, 2, 0], grid)
  c2 = CellField([0, 3, 4, 0], grid)
  c3 = CellField([0, 4, 6, 0], grid)
  c1_plus_c2 = c1 + c2
  c1_plus_c2.data == c3.data
end

function test_cell_times(T)
  nz = 2
  Lz = 2.0
  grid = UniformGrid(T, nz, Lz)
  c1 = CellField([0, 1, 2, 0], grid)
  c2 = CellField([0, 3, 4, 0], grid)
  c3 = CellField([0, 3, 8, 0], grid)
  c1_times_c2 = c1 * c2
  c1_times_c2.data == c3.data
end

function test_face_plus(T)
  nz = 2
  Lz = 2.0
  grid = UniformGrid(T, nz, Lz)
  c1 = FaceField([1, 2, 2], grid)
  c2 = FaceField([3, 4, 3], grid)
  c3 = FaceField([4, 6, 5], grid)
  c1_plus_c2 = c1 + c2
  c1_plus_c2.data == c3.data
end

function test_face_times(T)
  nz = 2
  Lz = 2.0
  grid = UniformGrid(T, nz, Lz)
  f1 = FaceField([1, 2, 2], grid)
  f2 = FaceField([3, 4, 3], grid)
  f3 = FaceField([3, 8, 6], grid)
  f1_times_f2 = f1 * f2
  f1_times_f2.data == f3.data
end

function test_diffusion_basic()
  model = Diffusion.Model(nz=4, Lz=2.0, κ=0.1)
  model.parameters.κ == 0.1
end

function test_diffusion_set_c()
  model = Diffusion.Model(nz=4, Lz=2.0, κ=0.1)
  c0 = [1.0, 2.0, 3.0, 4.0]
  Diffusion.set!(model, c=c0)
  model.solution.data[1:model.grid.nz] == c0
end

function test_diffusion()
  model = Diffusion.Model(nz=100, Lz=1.0, κ=0.1)
  model.parameters.κ == 0.1
end



# --
# Run tests
# --

@testset "Utils" begin
  @test test_zeros(Float64)
  @test test_zeros(Float32)
end

@testset "Grids" begin
  nz = 3
  Lz = 4.2
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
  nz = 3
  Lz = 4.2
  for T in (Float64, Float32, Float16)
    @test test_cell_field_construction(T, nz, Lz)
    @test test_face_field_construction(T, nz, Lz)
    @test test_cell_∂z(T)
    @test test_face_∂z(T)
    @test test_cell_plus(T)
    @test test_cell_times(T)
  end
  @test test_field_indexing()
end

@testset "Diffusion" begin
  @test test_diffusion_basic()
  @test test_diffusion_set_c()
end

#=
@testset "Forcing" begin

  function testforcinginterpolant(fld=:shortwave)
    kwargs = Dict(:tdata=>[0, year], fld=>[0, 1.0])
    forcing = Forcing(; kwargs...)
    finterp = ForcingInterpolant(forcing)
    isapprox(getfield(finterp, fld)(0.6year), 0.6)
  end

  for fld in fieldnames(ForcingInterpolant)
    @test testforcinginterpolant(fld)
  end
end
=#

#=
@testset "Mixing" begin
  include("testmixing.jl")
  @test testconvect1()
  @test testconvect2()
  @test testgradri()
  @test testbulkri()
  @test testbulkmix1()
  @test testbulkmix2()
  @test testgradmixing()
end

@testset "Time stepping" begin
  include("teststeps.jl")
  @test teststepuv1()
  @test teststepuv2()
  @test teststepforward()
end
=#
