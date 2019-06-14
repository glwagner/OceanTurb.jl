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
struct Equation{TM, TL, TK, TR, U}
    M       :: TM
    L       :: TL
    K       :: TK
    R       :: TR
    update! :: U
end

minuszero(args...) = -0
funkynothing(args...) = nothing

function Equation(;
         K,
         M = Tuple(minuszero for i=1:length(K)),
         L = Tuple(minuszero for i=1:length(K)),
         R = Tuple(minuszero for i=1:length(K)),
    update = funkynothing
    )
    return Equation(M, L, K, R, update)
end
