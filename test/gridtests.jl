# --
# Define tests
# --

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
