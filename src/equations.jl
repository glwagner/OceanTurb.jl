"""
    Equation(R, K)

Construct an Equation that represents a PDE of the form

    ∂ϕ/∂t = ∂/∂z (K ∂ϕ/∂z) + R ,

so that K is a diffusivity and R captures 'remaining' terms.

In OceanTurb, K and R are arrays of functions.
The functions K[j](m, i) and R[j](m, i) calculate the
diffusivity K and remaining terms R for the model `m` at
grid point `i`.
"""
struct Equation{TR, TK, U}
    R       :: TR
    K       :: TK
    update! :: U
end

nothing_func(args...) = nothing

Equation(R, K) = Equation(R, K, nothing_func)
