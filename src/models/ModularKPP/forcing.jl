addzero(args...) = 0

"""
    Forcing(; U=addzero, V=addzero, T=addzero, S=addzero)

Construct a `NamedTuple` of forcing functions for KPP `Model`s for each
field `U, V, T, S`. The functions must have the signature `forcing(model::Model, i)`,
where `i` is the vertical index at which the forcing is applied.
"""
Forcing(; U=addzero, V=addzero, T=addzero, S=addzero) = (U=U, V=V, T=T, S=S)

