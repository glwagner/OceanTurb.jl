using
  OceanTurb,
  Test

import OceanTurb

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

function test_uniform_grid_limits(T, nz, Lz)
  dz = Lz/nz
  grid = UniformGrid(T, nz, Lz)
  grid.zc[1] == -Lz+0.5*dz && grid.zf[1] == -Lz && grid.zc[nz] == -0.5*dz && grid.zf[nz] == 0
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
    #@test test_uniform_grid_limits(T, nz, Lz)
  end
end

@testset "Fields" begin
  nz = 3
  Lz = 4.2
  for T in (Float64, Float32, Float16)
    @test test_cell_field_construction(T, nz, Lz)
    @test test_face_field_construction(T, nz, Lz)
  end
  @test test_field_indexing()
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
