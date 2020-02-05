#
# Define tests
# 

function fill_ghost_cells!(c, κtop, κbottom, model, fieldbcs)
    fill_bottom_ghost_cell!(fieldbcs.bottom, c, κbottom, model)
    fill_top_ghost_cell!(fieldbcs.top, c, κtop, model, c.grid.N)
    return nothing
end

function cell_field_construction(T, N, L)
    grid = UniformGrid(T, N, L)
    c = CellField(grid)
    typeof(c) <: CellField
end

function face_field_construction(T, N, L)
    grid = UniformGrid(T, N, L)
    f = FaceField(grid)
    typeof(f) <: FaceField
end

function field_indexing()
    grid = UniformGrid(3, 4.2)
    c = CellField(grid)
    f = FaceField(grid)
    val = 2.1
    c[2] = val
    f[2] = val
    c[2] == val && f[2] == val
end

function cell_∂z(T)
    grid = UniformGrid(T, 2, 2.0)
    c = CellField([2, 4], grid)
    cz = ∂z(c)
    cz_answer = FaceField([0, 2, 0], grid)
    cz.data == cz_answer.data
end

function face_∂z(T)
    grid = UniformGrid(T, 2, 2.0)
    f = FaceField([2, 4, 6], grid)
    fz = ∂z(f)
    fz_answer = CellField([2, 2], grid)
    fz.data[1:2] == fz_answer.data[1:2]
end

function cell_plus(T)
    grid = UniformGrid(T, 4, 2.0)
    c1 = CellField([0, 1, 2, 0], grid)
    c2 = CellField([0, 3, 4, 0], grid)
    c3 = CellField([0, 4, 6, 0], grid)
    c1_plus_c2 = c1 + c2
    c1_plus_c2.data == c3.data
end

function cell_times(T)
    grid = UniformGrid(T, 4, 2.0)
    c1 = CellField([0, 1, 2, 0], grid)
    c2 = CellField([0, 3, 4, 0], grid)
    c3 = CellField([0, 3, 8, 0], grid)
    c1_times_c2 = c1 * c2
    c1_times_c2.data == c3.data
end

function face_plus(T)
    grid = UniformGrid(T, 2, 2.0)
    c1 = FaceField([1, 2, 2], grid)
    c2 = FaceField([3, 4, 3], grid)
    c3 = FaceField([4, 6, 5], grid)
    c1_plus_c2 = c1 + c2
    c1_plus_c2.data == c3.data
end

function face_times(T)
    grid = UniformGrid(T, 2, 2.0)
    f1 = FaceField([1, 2, 2], grid)
    f2 = FaceField([3, 4, 3], grid)
    f3 = FaceField([3, 8, 6], grid)
    f1_times_f2 = f1 * f2
    f1_times_f2.data == f3.data
end

function set_scalar_field(loc, T)
    g = UniformGrid(T, 2, 2.0)
    f = Field(loc, g)
    a = 7
    set!(f, a)
    !any(@. !(f.data == a))
end

function set_array_field(loc, T)
    g = UniformGrid(T, 2, 2.0)
    f = Field(loc, g)
    data = rand(1:10, length(f))
    set!(f, data)
    f.data[1:length(f)] == data
end

function set_function_field(loc, T)
    g = UniformGrid(T, 2, 2.0)
    f = Field(loc, g)
    fcn(z) = z^2
    z = nodes(f)
    set!(f, fcn)
    i = 2
    f[i] == fcn(z[i])
end

function integral_is_correct(T)
    grid = UniformGrid(T, 3, 3.0)
    c = CellField([1, 2, 3], grid)
    integral(c) / grid.L == 2
end

function integral_range(T, N, L, z₋, z₊)
    grid = UniformGrid(T, N, L)
    c = CellField(ones(N), grid)
    integral(c, z₋, z₊) == z₊ - z₋
end


struct DummyModel end

function ghost_cell_value(T)
    model = DummyModel()
    grid = UniformGrid(T, 3, 3.0)
    c = CellField([1, 2, 3], grid)
    c_top = T(0.2)
    c_bottom = T(0.11)
    top_bc = ValueBoundaryCondition(c_top)
    bottom_bc = ValueBoundaryCondition(c_bottom)
    fieldbcs = FieldBoundaryConditions(bottom_bc, top_bc)
    fill_ghost_cells!(c, 0, 0, model, fieldbcs)

    onface(c, 1) ≈ c_bottom && onface(c, grid.N+1) ≈ c_top
end

function ghost_cell_gradient(T)
    model = DummyModel()
    grid = UniformGrid(T, 3, 3.0)
    c = CellField([1, 2, 3], grid)
    cz_top = T(0.2)
    cz_bottom = T(0.11)
    top_bc = GradientBoundaryCondition(cz_top)
    bottom_bc = GradientBoundaryCondition(cz_bottom)
    fieldbcs = FieldBoundaryConditions(bottom_bc, top_bc)

    fill_ghost_cells!(c, 0, 0, model, fieldbcs)

    ∂z(c, 1) ≈ cz_bottom && ∂z(c, grid.N+1) ≈ cz_top
end

function ghost_cell_flux(T)
    model = DummyModel()
    grid = UniformGrid(T, 3, 3.0)
    c = CellField([1, 2, 3], grid)
    Fc_top = T(0.2)
    Fc_bottom = T(0.11)
    top_bc = FluxBoundaryCondition(Fc_top)
    bottom_bc = FluxBoundaryCondition(Fc_bottom)
    fieldbcs = FieldBoundaryConditions(bottom_bc, top_bc)

    κ = T(0.47)
    fill_ghost_cells!(c, κ, κ, model, fieldbcs)

    (OceanTurb.flux(0, κ, c, 1) ≈ Fc_bottom
        && OceanTurb.flux(0, κ, c, grid.N+1) ≈ Fc_top)
end

function absolute_error_is_correct(T)
    grid = UniformGrid(T, 3, 3)
    c = CellField([0, 0, 0], grid)
    d = CellField([2, 2, 2], grid)
    absolute_error(c, d) == 2
end

function absolute_error_fine_to_coarse()
    fine_grid = UniformGrid(N=4, L=1)
    coarse_grid = UniformGrid(N=2, L=1)
    fine_field = CellField([2, 4, 6, 8], fine_grid)
    coarse_field = CellField([0, 1], coarse_grid)
    return absolute_error(coarse_field, fine_field) == sqrt(45/2)
end

function absolute_error_coarse_to_fine()
    fine_grid = UniformGrid(N=4, L=1)
    coarse_grid = UniformGrid(N=2, L=1)
    fine_field = CellField([2, 4, 6, 8], fine_grid)
    coarse_field = CellField([0, 1], coarse_grid)
    return absolute_error(fine_field, coarse_field) == sqrt(47/2)
end

function relative_error_is_correct(T)
    grid = UniformGrid(T, 3, 1)
    c = CellField([0, 0, 0], grid)
    d = CellField([2, 2, 2], grid)
    relative_error(c, d) == 1
end

@testset "Fields" begin
    nz, Lz = 3, 4.2
    for T in (Float64,)
        @test cell_field_construction(T, nz, Lz)
        @test face_field_construction(T, nz, Lz)
        @test cell_∂z(T)
        @test face_∂z(T)
        @test cell_plus(T)
        @test cell_times(T)
        for loc in (Face, Cell)
            @test set_scalar_field(loc, T)
            @test set_array_field(loc, T)
            @test set_function_field(loc, T)
        end
        @test integral_is_correct(T)
        @test integral_range(T, 10, 10, -5, -3)
        @test integral_range(T, 10, 10, -5.1, -4.5)
        @test integral_range(T, 17, 3, -2.5, -1.3)
        @test integral_range(T, 3, 30, -2.5, -1.3)
        @test ghost_cell_value(T)
        @test ghost_cell_gradient(T)
        @test ghost_cell_flux(T)
        @test absolute_error_is_correct(T)
        @test relative_error_is_correct(T)
    end
    @test absolute_error_fine_to_coarse()
    @test absolute_error_coarse_to_fine()
    @test field_indexing()
end

