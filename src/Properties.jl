module Properties

  export ConformalMap, DerivativeMap

  abstract type ConformalMap end

  struct DerivativeMap{M<:ConformalMap}
    m :: M
  end

end
