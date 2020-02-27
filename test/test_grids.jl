function test_uniform_grid_type(T, N, H)
  grid = UniformGrid(T, N, H)
  T == typeof(grid.H) && T == eltype(grid.zc) && T == eltype(grid.zf)
end

function test_uniform_grid_spacing(T, N, H)
  Δ = T(H/N)
  grid = UniformGrid(T, N, H)
  Δ == grid.Δc && Δ == grid.Δf
end

function test_uniform_grid_limits_zc_left(T, N, H)
  Δ = H/N
  grid = UniformGrid(T, N, H)
  grid.zc[1] ≈ -H+0.5Δ
end

function test_uniform_grid_limits_zc_right(T, N, H)
  Δ = H/N
  grid = UniformGrid(T, N, H)
  grid.zc[N] ≈ -0.5Δ
end

function test_uniform_grid_limits_zf_left(T, N, H)
  Δ = H/N
  grid = UniformGrid(T, N, H)
  grid.zf[1] ≈ -H
end

function test_uniform_grid_limits_zf_right(T, N, H)
  Δ = H/N
  grid = UniformGrid(T, N, H)
  grid.zf[N+1] ≈ 0
end

@testset "Grids" begin
    N, H = 3, 4.2
    for T in (Float64, Float32, Float16)
        @test test_uniform_grid_type(T, N, H)
        @test test_uniform_grid_spacing(T, N, H)
        @test test_uniform_grid_limits_zc_right(T, N, H)
        @test test_uniform_grid_limits_zf_right(T, N, H)
        @test test_uniform_grid_limits_zc_left(T, N, H)
        @test test_uniform_grid_limits_zf_left(T, N, H)
    end
end

