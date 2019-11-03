struct StandardCubicPolynomial <: AbstractParameters end

Base.@kwdef struct GeneralizedCubicPolynomial{T} <: AbstractParameters
    CS0 :: T = 0.0
    CS1 :: T = 1.0
end

shape(d, p::StandardCubicPolynomial) = d * (1-d)
shape(d, p::GeneralizedCubicPolynomial) = d * (1-d) * ( p.CS0 + p.CS1*(1-d) )
