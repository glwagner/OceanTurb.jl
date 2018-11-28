function testzeros(T=Float64, dims=(13, 45))
  a1, b1 = zeros(T, dims), zeros(T, dims)
  PWP.@zeros T dims a2 b2
  eltype(a1) == T && a1 == a2 && b1 == b2
end

function testexample(H=200, nz=200)
  examplemodel = loadexample(H, nz)
  nz == examplemodel.ocean.nz
end

function testparams()
  params = Parameters(f=π)
  isapprox(params.f, π)
end
