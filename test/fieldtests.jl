# --
# Define tests
# --

function test_cell_field_construction(T, nz, Lz)
  grid = UniformGrid(T, nz, Lz)
  c = CellField(grid)
  length(c.data) == OceanTurb.cell_length(nz) && size(c.data) == OceanTurb.cell_size(nz)
end

function test_face_field_construction(T, nz, Lz)
  grid = UniformGrid(T, nz, Lz)
  f = FaceField(grid)
  length(f.data) == OceanTurb.face_length(nz) && size(f.data) == OceanTurb.face_size(nz)
end

function test_field_indexing()
  grid = UniformGrid(3, 4.2)
  c = CellField(grid)
  f = FaceField(grid)
  val = 2.1
  c[2] = val
  f[2] = val
  c[2] == val && f[2] == val
end

function test_cell_∂z(T)
  grid = UniformGrid(T, 2, 2.0)
  c = CellField([2, 4], grid)
  cz = ∂z(c)
  cz_answer = FaceField([0, 2, 0], grid)
  cz.data == cz_answer.data
end

function test_face_∂z(T)
  grid = UniformGrid(T, 2, 2.0)
  f = FaceField([2, 4, 6], grid)
  fz = ∂z(f)
  fz_answer = CellField([2, 2], grid)
  fz.data == fz_answer.data
end

function test_cell_plus(T)
  grid = UniformGrid(T, 4, 2.0)
  c1 = CellField([0, 1, 2, 0], grid)
  c2 = CellField([0, 3, 4, 0], grid)
  c3 = CellField([0, 4, 6, 0], grid)
  c1_plus_c2 = c1 + c2
  c1_plus_c2.data == c3.data
end

function test_cell_times(T)
  grid = UniformGrid(T, 4, 2.0)
  c1 = CellField([0, 1, 2, 0], grid)
  c2 = CellField([0, 3, 4, 0], grid)
  c3 = CellField([0, 3, 8, 0], grid)
  c1_times_c2 = c1 * c2
  c1_times_c2.data == c3.data
end

function test_face_plus(T)
  grid = UniformGrid(T, 2, 2.0)
  c1 = FaceField([1, 2, 2], grid)
  c2 = FaceField([3, 4, 3], grid)
  c3 = FaceField([4, 6, 5], grid)
  c1_plus_c2 = c1 + c2
  c1_plus_c2.data == c3.data
end

function test_face_times(T)
  grid = UniformGrid(T, 2, 2.0)
  f1 = FaceField([1, 2, 2], grid)
  f2 = FaceField([3, 4, 3], grid)
  f3 = FaceField([3, 8, 6], grid)
  f1_times_f2 = f1 * f2
  f1_times_f2.data == f3.data
end

function test_set_scalar_field(loc, T)
  g = UniformGrid(T, 2, 2.0)
  f = Field(loc, g)
  a = 7
  set!(f, a)
  f.data == a*ones(length(f))
end

function test_set_array_field(loc, T)
  grid = UniformGrid(T, 2, 2.0)
  f = Field(loc, grid)
  data = rand(1:10, length(f))
  set!(f, data)
  f.data == data
end

function test_set_function_field(loc, T)
  grid = UniformGrid(T, 2, 2.0)
  f = Field(loc, grid)
  fcn(z) = z^2
  data_ans = fcn.(zdata(f))
  set!(f, fcn)
  f.data == data_ans
end
