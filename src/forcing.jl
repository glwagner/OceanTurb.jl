alternatenametable = Dict(
  "τˣ"=>"stress_x",
  "τʸ"=>"stress_y"
)

struct Forcing{T}
  ndata::Int
  tdata::AbstractArray{T,1}
  shortwave::AbstractArray{T,1}
  longwave::AbstractArray{T,1}
  latentflux::AbstractArray{T,1}
  sensibleflux::AbstractArray{T,1}
  τˣ::AbstractArray{T,1}
  τʸ::AbstractArray{T,1}
  precip::AbstractArray{T,1}
  evap::AbstractArray{T,1}
end

function Forcing(; tdata=[0, year], 
                  shortwave=0tdata, longwave=0tdata, latentflux=0tdata, 
                  sensibleflux=0tdata, stress_x=0tdata, stress_y=0tdata, 
                  precip=0tdata, evap=0tdata,
                  τˣ=stress_x, τʸ=stress_y
                 )

  Forcing(length(tdata), tdata, shortwave, longwave, latentflux, sensibleflux, 
          τˣ, τʸ, precip, evap)
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
    forcingfields["latentflux"],
    forcingfields["sensibleflux"],
    forcingfields["τˣ"],
    forcingfields["τʸ"],
    forcingfields["precip"],
    forcingfields["evap"]
   )
end

struct ForcingInterpolant{T}
  shortwave::ScaledInterpolation{T}
  longwave::ScaledInterpolation{T}
  latentflux::ScaledInterpolation{T}
  sensibleflux::ScaledInterpolation{T}
  τˣ::ScaledInterpolation{T}
  τʸ::ScaledInterpolation{T}
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
