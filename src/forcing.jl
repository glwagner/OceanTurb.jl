struct Forcing{T}
  ndata::Int
  tdata::AbstractArray{T,1}
  shortwave::AbstractArray{T,1}
  longwave::AbstractArray{T,1}
  latent::AbstractArray{T,1}
  sensible::AbstractArray{T,1}
  xstress::AbstractArray{T,1}
  ystress::AbstractArray{T,1}
  precip::AbstractArray{T,1}
  evap::AbstractArray{T,1}
end

"""
    Forcing(; tdata=[0, year], forcingfields...)

Construct a `Forcing` specified at time points in `tdata`. The forcing
fields, which are arrays of data corresponding to the times in `tdata` and 
are assumed to be measured/calculated/specified at the surface
and `z=0`, are specified by keyword arguments.
Their default values are `0tdata`. Possible forcing inputs are

  * `shortwave` : incoming shortwave radiation
  * `longwave`  : outgoing longwave radiation
  * `latent`    : incoming latent heat flux
  * `sensible`  : incoming latent heat flux
  * `precip`    : precipitation
  * `evap`      : evaporation
  * `xstress`   : wind stress in the `x`-direction
  * `ystress`   : wind stress in the `y`-direction

"""
function Forcing(; tdata=[0, year], 
                  shortwave=0tdata, longwave=0tdata,
                  latent=0tdata, sensible=0tdata,
                  precip=0tdata, evap=0tdata,
                  xstress=0tdata, ystress=0tdata
                 )

  Forcing(length(tdata), tdata, shortwave, longwave, latent, sensible, 
          xstress, ystress, precip, evap)
end

iskey(key, c) = key in keys(c)
setforcingvar(varsym, forcingdict) = eval(varsym) = forcingdict[String(varsym)]

"""
    loadforcing(filepath)

Initalize a `Forcing` from the JLD2 file at `filepath`. The names of the file's 
fields must correspond to the names of the `Forcing` fields.
"""
function loadforcing(filepath)

  forcingfields = Dict{String,Any}()

  file = jldopen(filepath, "r")
  try
    iskey("t", file) || error("The JLD2 file $filepath has no `t` key!")

    t = file["t"]
    nt = length(t)
    forcingfields["nt"] = nt
    forcingfields["t"] = t

    # Extract rest of fields
    for fieldsym in fieldnames(Forcing)
      if fieldsym != :nt
        fieldstr = String(fieldsym)
        filefieldstr = (!iskey(fieldstr, file) && 
                        iskey(fieldstr, alternatenametable)
                        ? alternatenametable[fieldstr] : fieldstr)
        forcingfields[fieldstr] = (iskey(filefieldstr, file) ? 
                                   file[filefieldstr] : 0t)
        nt == length(forcingfields[fieldstr]) || (
          error("Length of forcing field $fieldstr is incompatable with t"))
      end
    end

  finally
    close(file)
  end

  Forcing(
    forcingfields["nt"],
    forcingfields["t"],
    forcingfields["shortwave"],
    forcingfields["longwave"],
    forcingfields["latent"],
    forcingfields["sensible"],
    forcingfields["xstress"],
    forcingfields["ystress"],
    forcingfields["precip"],
    forcingfields["evap"]
   )
end

struct ForcingInterpolant{T}
  shortwave::ScaledInterpolation{T}
  longwave::ScaledInterpolation{T}
  latent::ScaledInterpolation{T}
  sensible::ScaledInterpolation{T}
  xstress::ScaledInterpolation{T}
  ystress::ScaledInterpolation{T}
  precip::ScaledInterpolation{T}
  evap::ScaledInterpolation{T}
end

function ForcingInterpolant(forcing)
  dt = forcing.tdata[2]-forcing.tdata[1]
  t⁰ = forcing.tdata[1]
  tᶠ = forcing.tdata[end]
  for name in fieldnames(ForcingInterpolant)
    field = getfield(forcing, name)
    @eval begin
      $name = scale(interpolate($field, BSpline(Linear())), $t⁰:$dt:$tᶠ)
    end
  end
  ForcingInterpolant(eval.(fieldnames(ForcingInterpolant))...)
end
