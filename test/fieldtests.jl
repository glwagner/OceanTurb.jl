# --
# Define tests
# --

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

