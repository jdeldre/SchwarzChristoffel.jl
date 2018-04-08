module Properties

  export ConformalMap, DerivativeMap

  abstract type ConformalMap end

  struct DerivativeMap{M<:ConformalMap}
    m :: M
  end
  function Base.show(io::IO, dm::DerivativeMap)
      println(io, "d/dζ of $(dm.m)")
  end

end
