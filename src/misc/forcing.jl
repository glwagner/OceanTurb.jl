iskey(key, c) = key in keys(c)
setforcingvar(varsym, forcingdict) = eval(varsym) = forcingdict[String(varsym)]

struct ForcingData{A}
  ndata::Int
  tdata::A
  shortwave::A
  longwave::A
  latent::A
  sensible::A
  xstress::A
  ystress::A
  precip::A
  evap::A
end

struct ForcingInterpolant{T}
  shortwave::ScaledInterpolation{T}
  longwave::ScaledInterpolation{T}
  latent::ScaledInterpolation{T}
  sensible::ScaledInterpolation{T}
  precip::ScaledInterpolation{T}
  evap::ScaledInterpolation{T}
  xstress::ScaledInterpolation{T}
  ystress::ScaledInterpolation{T}
end

struct Forcing{T,A}
  data::ForcingData{A}
  interp::ForcingInterpolant{T}
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
function ForcingData(; 
      tdata = [0, year], 
  shortwave = 0tdata, 
   longwave = 0tdata,
     latent = 0tdata, 
   sensible = 0tdata,
     precip = 0tdata, 
       evap = 0tdata,
    xstress = 0tdata, 
    ystress = 0tdata
)

  dt = tdata[2]-tdata[1]
  if any(dt .!= tdata[2:end]-tdata[1:end-1]) 
    error("Data must be evenly spaced in time")
  end

  ndata = length(tdata) 
  # TODO: ensure all fields have the same length
  # for name in fieldnames(ForcingInterpolant)...
  
  ForcingData(tdata, shortwave, longwave, latent, sensible, xstress, ystress, 
              precip, evap)
              
end

function ForcingInterpolant(forcing)
  dt = forcing.tdata[2]-forcing.tdata[1]
  t0 = forcing.tdata[1]
  tf = forcing.tdata[end]
  for name in fieldnames(ForcingInterpolant)
    field = getfield(forcing, name)
    @eval begin
      $name = scale(interpolate($field, BSpline(Linear())), $t0:$dt:$tf)
    end
  end
  ForcingInterpolant(eval.(fieldnames(ForcingInterpolant))...)
end

#=
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
=#
