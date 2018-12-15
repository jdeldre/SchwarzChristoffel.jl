module MapTypes

export ConformalMap, DerivativeMap, InverseMap

abstract type ConformalMap end

struct DerivativeMap{M<:ConformalMap}
  m :: M
end
function Base.show(io::IO, dm::DerivativeMap)
  println(io, "d/dζ of $(dm.m)")
end

struct InverseMap{M<:ConformalMap}
  m :: M
end
function Base.show(io::IO, minv::InverseMap)
  println(io, "Inverse of $(minv.m)")
end

end
