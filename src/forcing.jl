struct Forcing{T}
  nt::Int
  t::T
  shortwave::T
  longwave::T
  latentflux::T
  sensibleflux::T
  stress_x::T
  stress_y::T
  precip::T
  evap::T
end

function Forcing(t; shortwave=nothing, longwave=nothing, latentflux=nothing, sensibleflux=nothing,
                    stress_x=nothing, stress_y=nothing, precip=nothing, evap=nothing)

  for fld in fieldnames(Forcing)
    eval(fld) == nothing && @eval $fld = 0t
    length(t) == length(eval(fld)) || error("All forcing fields must have length(t)!")
  end

  Forcing(length(t), eval.(fieldnames(Forcing)))
end

iskey(key, c) = key in keys(c)

setforcingvar(varsym, forcingdict) = eval(varsym) = forcingdict[String(varsym)]

"""
    loadforcing(filepath)

Initalize a `Forcing` from the JLD2 file at `filepath`. The names of the file's fields must
correspond to the names of the `Forcing` fields.
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
        forcingfields[fieldstr] = iskey(fieldstr, file) ? file[fieldstr] : 0t
        nt == length(forcingfields[fieldstr]) || error("Length of forcing field $fieldstr is incompatable with t")
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
    forcingfields["stress_x"],
    forcingfields["stress_y"],
    forcingfields["precip"],
    forcingfields["evap"]
   )
end
