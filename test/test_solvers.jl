function test_tridiagonalsolve()
    # Test problem:
    # A = [[1 4 0 0]; [3 4 1 0]; [0 2 3 4]; [0 0 1 3]]
    # x = [1; 2; 3; 4]
    # = > A*x = d = [9; 14; 29; 15]

    a = [3.0, 2.0, 1.0]
    b = [1.0, 4.0, 3.0, 3.0]
    c = [4.0, 1.0, 4.0]
    d = [9.0, 14.0, 29.0, 15.0]

    x = zeros(length(b))
    OceanTurb.tridiagonalsolve!(x, d, a, b, c)

    solution = [1.0, 2.0, 3.0, 4.0]
    x == solution
end

function test_ldiv(N)
    x = rand(N)
    A = Tridiagonal(rand(N, N))
    b = A*x
    test_x = zeros(N)

    ldiv!(test_x, A, b)

    x ≈ test_x
end

@testset "Solvers" begin
    @test test_tridiagonalsolve()
    @test test_ldiv(10)
end

