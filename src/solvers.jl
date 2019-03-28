import LinearAlgebra: ldiv!

"""
    tridiagonalsolve!(x, b, u, d, l)

Solve the tridiagonal system `A*x = b`, where `A` is tridiagonal and defined
by the upper diagonal `u`, diagonal `d`, and lower diagonal `l`

Following https://julialang.org/blog/2016/02/iteration for the multidimensional implementation.
"""
function ldiv!(x, A::Tridiagonal, b)
    N = length(x)
    N == length(b) || throw("x and b must have the same length.")
    _tridiagonalsolve!(x, b, A.dl, A.d, A.du, N)
    return nothing
end

function tridiagonalsolve!(x, d, a, b, c)
    length(x) == length(d) == (length(a)+1) == length(b) == (length(c)+1) || (
        throw(ArgumentError("The size of x, d, a, b, and c do not match.")))
    _tridiagonalsolve!(x, d, a, b, c, length(x))
end

function _tridiagonalsolve!(x, d, a, b, c, N)
    for i = 2:N
        @inbounds begin
            b[i] -= a[i-1] * c[i-1] / b[i-1]
            d[i] -= a[i-1] * d[i-1] / b[i-1]
        end
    end

    @inbounds x[N] = d[N] / b[N]

    for i = N-1:-1:1
        # Solve for x by bulk substitution
        @inbounds x[i] = (d[i] - c[i]*x[i+1] ) / b[i]
    end

    return nothing
end
