"""
    ModelBoundaryConditions([FT=Float64;] U = DefaultBoundaryConditions(FT),
                                          V = DefaultBoundaryConditions(FT),
                                          T = DefaultBoundaryConditions(FT),
                                          S = DefaultBoundaryConditions(FT))

Returns a `NamedTuple` of boundary conditions for a `KPP.Model` with solution
fields `U`, `V`, `T`, `S`.

Example
=======

julia> surface_temperature_flux(model) = cos(model.clock.time)

julia> T_bcs = BoundaryConditions(top = FluxBoundaryCondition(surface_flux))

julis> bcs = ModularKPP.ModelBoundaryConditions(T=T_bcs)
"""
function ModelBoundaryConditions(FT=Float64; U = DefaultBoundaryConditions(FT),
                                             V = DefaultBoundaryConditions(FT),
                                             T = DefaultBoundaryConditions(FT),
                                             S = DefaultBoundaryConditions(FT))
    return (U=U, V=V, T=T, S=S)
end

