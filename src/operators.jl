#
# Local diffusive flux
#

const BC = BoundaryCondition

# ∇K∇c for c::CellField
K∂z(K, c, i) = K*∂z(c, i)
∇K∇c(Kᵢ₊₁, Kᵢ, c, i)            = ( K∂z(Kᵢ₊₁, c, i+1) -    K∂z(Kᵢ, c, i)      ) /    Δf(c, i)
∇K∇c_top(Kᴺ, c, top_flux)       = (     -top_flux     - K∂z(Kᴺ, c, length(c)) ) / Δf(c, length(c))
∇K∇c_bottom(K₂, c, bottom_flux) = (   K∂z(K₂, c, 2)  +      bottom_flux       ) /    Δf(c, 1)

## Top and bottom flux estimates for constant (Dirichlet) boundary conditions
bottom_flux(K, c, c_bndry, Δf) = -2K*( bottom(c) - c_bndry ) / Δf # -K*∂c/∂z at the bottom
top_flux(K, c, c_bndry, Δf)    = -2K*(  c_bndry  -  top(c) ) / Δf # -K*∂c/∂z at the top

∇K∇c_top(Kᴺ⁺¹, Kᴺ, c, bc::BC{<:Flux}, model) = ∇K∇c_top(Kᴺ, c, getbc(model, bc))
∇K∇c_bottom(K₂, K₁, c, bc::BC{<:Flux}, model) = ∇K∇c_bottom(K₂, c, getbc(model, bc))
∇K∇c_bottom(K₂, K₁, c, bc::BC{<:Gradient}, model) = ∇K∇c_bottom(K₂, c, -K₁*getbc(model, bc))

function ∇K∇c_bottom(K₂, K₁, c, bc::BC{<:Value}, model)
    flux = bottom_flux(K₁, c, getbc(model, bc), Δf(model.grid, 1))
    return ∇K∇c_bottom(K₂, c, flux)
end

function ∇K∇c_top(Kᴺ⁺¹, Kᴺ, c, bc::BC{<:Value}, model)
    flux = top_flux(Kᴺ⁺¹, c, getbc(model, bc), Δf(model.grid, length(c)+1))
    return ∇K∇c_top(Kᴺ, c, flux)
end
