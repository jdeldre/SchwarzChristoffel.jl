module Exterior

using NLsolve

using ..Properties
using ..Polygons
using ..Integration

include("Reindex.jl")
using .Reindex

export PowerMap,ExteriorMap,parameters,coefficients,
        moments,area,centroid,Jmoment


struct PowerMap <: ConformalMap
    "power series coefficients, ccoeff[1] -> c₁, ccoeff[2] -> c₀, ccoeff[3] -> c₋₁, etc"
    ccoeff::Vector{Complex128}

    "number of multipole coefficients (in addition to c₁ and c₀)"
    ncoeff::Int

    "number of plotting control points"
    N::Int

    "control point coordinates in circle space"
    ζ::Vector{Complex128}

    "control point coordinates in body-fixed space"
    z::Vector{Complex128}

    "map Jacobian in body-fixed coordinates"
    dzdζ::Vector{Complex128}

    "coefficients of power series of |z̃(ζ)|²"
    dcoeff::Vector{Complex128}

    "Area enclosed by the mapped shape"
    area      :: Float64

    "Centroid of the mapped shape"
    Zc        :: Complex128

    "2nd area moment of the mapped shape"
    J         :: Float64

end

circle(N) = [exp(im*2π*(i-1)/N) for i in 1:N]

"""
    PowerMap(c::Vector{Complex12}[;N = 200]) <: ConformalMap

Create a power series map from the exterior of the unit
circle to the exterior of a shape defined by the power series coefficients.

# Example

```jldoctest
julia> c = Complex128[1,0,1/4];

julia> m = PowerMap(c)
Power series map:
   multipole coefficients: c₁ = 1.0 + 0.0im, c₀ = 0.0 + 0.0im, c₋ᵢ = 0.25 + 0.0im, i = 1:1
```
"""
function PowerMap(ccoeff::Vector{Complex128}; N::Int = 200)
  # Must have at least two entries. If no entries, make leading entry
  # a 1 (an identity map)
  if isempty(ccoeff)
    push!(ccoeff,1)
  end
  # if only one entry, add a zero for c₀
  if length(ccoeff) == 1
    push!(ccoeff,0)
  end
  ncoeff = length(ccoeff)-2

  ζ = circle(N)
  z = powerseries(ζ,ccoeff)
  dz, ddz = d_powerseries(ζ,ccoeff)


  # Coefficients of |z(ζ)|²
  # dcoeff[1] = d₀, dcoeff[2] = d₋₁, etc.
  # Note that d₁ = conj(d₋₁) = conj(dcoeff[2]), d₂ = conj(d₋₂) = conj(dcoeff[3])
  dcoeff = [dot(ccoeff,ccoeff)]
  for k = 1:ncoeff+1
      push!(dcoeff,dot(ccoeff[1:end-k],ccoeff[k+1:end]))
  end

  c = Reflect(ShiftReindex(ccoeff,-2))
  d = Reflect(OddReindex(dcoeff))

  k = -1:ncoeff
  l = -1:ncoeff
  kml = k[:,ones(Int,ncoeff+2)]'-l[:,ones(Int,ncoeff+2)]
  area = -π*sum(k.*abs.(c(-k)).^2)

  if area > 0
    Zc = -π/area*sum(k.*c(-k).*d(k))
  else
    Zc = mean(z)
  end
  J = Float64(-0.5π*c(-k)'*d(-kml)*(l.*c(-l)))


  PowerMap(ccoeff, ncoeff, N, ζ, z, dz, dcoeff, area, Zc, J)
end

function Base.show(io::IO, m::PowerMap)
    println(io, "Power series map:")
    print(io, "   multipole coefficients: c₁ = $(m.ccoeff[1]), ")
    print(io, "c₀ = $(m.ccoeff[2]), ")
    print(io,"c₋ᵢ = ")
    for i = 1:m.ncoeff
      print(io,"$(m.ccoeff[2+i]), ")
    end
    println(io, "i = 1:$(m.ncoeff)")
end

"""
    m::PowerMap(ζ::T) -> T

Evaluates the power-series mapping `m` at the points `ζ`.

# Example

```jldoctest
julia> c = Complex128[1,0,1/4];

julia> m = PowerMap(c);

julia> ζ = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];

julia> m(ζ)
3-element Array{Complex{Float64},1}:
   0.595+1.515im
 -1.2125-0.9875im
     0.0+0.195909im
```
"""
(m::PowerMap)(ζ) = powerseries(ζ,m.ccoeff)

evalderiv(ζ,m::PowerMap) = d_powerseries(ζ,m.ccoeff)

# struct PowerSeries{T}
#   C :: Vector{T}
# end
# function (f::PowerSeries{T})(ζ::Number) where T
#   ζⁿ = ζ
#   z = zero(ζ)
#   for c in f.C
#     z += c*ζⁿ
#     ζⁿ /= ζ
#   end
#   return z
# end
#
# (f::PowerSeries{T})(ζ::Vector{Number}) = f.(ζ)

function powerseries(ζ::Complex128,C::Vector{Complex128})
  ζⁿ = ζ
  z = zero(ζ)
  for c in C
    z += c*ζⁿ
    ζⁿ /= ζ
  end
  z
end

powerseries(ζs::Vector{Complex128},C::Vector{Complex128}) =
              [powerseries(ζ,C) for ζ in ζs]

function d_powerseries(ζ::Complex128,C::Vector{Complex128})
  dz = C[1]
  ζⁿ = 1/ζ^2
  ddz = Complex128(0)
  for n in 1:length(C)-2
    dz -= n*C[n+2]*ζⁿ
    ζⁿ /= ζ
    ddz += n*(n+1)*ζⁿ
  end
  return dz, ddz
end

function d_powerseries(ζs::Vector{Complex128},C::Vector{Complex128})

  dz = zeros(ζs)
  ddz = zeros(ζs)
  for (i,ζ) in enumerate(ζs)
    dz[i], ddz[i] = d_powerseries(ζ,C)
  end
  return dz, ddz

end




#=   Exterior map from to polygon  =#

struct ExteriorMap <: ConformalMap

  "Number of vertices on polygon"
  N :: Int

  "Coordinates of vertices on the polygon, defined ccw"
  z :: Vector{Complex128}

  "Interior angles in the polygon"
  angle  :: Vector{Float64}

  "Quadrature nodes and weights for the Gauss-Jacobi"
  qdata  :: Tuple{Array{Float64,2},Array{Float64,2}}

  "Accuracy of the quadrature"
  accuracy  :: Float64

  "Coordinates of the prevertices on the unit circle (interior)"
  ζ :: Vector{Complex128}

  "Constant factor of the mapping"
  constant  :: Complex128

  "Coordinates of the pre-prevertices on the unit circle (exterior)"
  preprev   :: Vector{Complex128}

  "Angular locations of the pre-prevertices on the unit circle (exterior)"
  prevangle :: Vector{Float64}

  "Number of multipole coefficients"
  ncoeff    :: Int64

  "Multipole coefficients for z(ζ). ccoeff[1] -> c₁, ccoeff[2] -> c₀, etc."
  ccoeff    :: Vector{Complex128}

  "Multipole coefficients for |z(ζ)|². dcoeff[1] -> d₀, doceff[2] -> d₋₁, etc."
  dcoeff    :: Vector{Complex128}

  "Moments of prevertices"
  mom       :: Vector{Complex128}

  "Area enclosed by the mapped polygon"
  area      :: Float64

  "Centroid of the mapped polygon"
  Zc        :: Complex128

  "2nd area moment of the mapped polygon"
  J         :: Float64

end

"""
    ExteriorMap(p::Polygon[;tol::Float64][,ncoeff::Int]) <: ConformalMap

Create a Schwarz-Christoffel map from the interior or exterior of
the unit circle to the exterior of polygon `p`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p)
Exterior map with 4 vertices
   vertices: (-1.0,-1.0), (0.2,-1.0), (1.0,0.5), (-1.0,1.0),
   interior angles/π: 0.5, 0.656, 0.422, 0.422,
   prevertices on circle: (1.0,0.0), (0.3764,-0.9265), (-0.9024,-0.4309), (-0.1868,0.9824),
   prevertex angles/π: -0.7291, -0.3519, 0.1291, 0.7111,
   constant = 0.6722 + 0.7669im, accuracy = 1.0e-8,
   number of multipole coefficients = 12
```

`ExteriorMap(p;tol=1e-12)` manually sets the tolerance to `1e-12`
(the default is 1e-8).

`ExteriorMap(p;ncoeff=16)` manually sets the number of coefficients of
negative powers of the multipole expansion of the mapping to `16`
(the default is 12).
"""
function ExteriorMap(p::Polygon;tol::Float64 = 1e-8,ncoeff::Int = 12)

  n = length(p.vert)

  zeta0 = Complex128[]

  w = flipdim(vertex(p),1)
  beta = 1.-flipdim(interiorangle(p),1)

  # do some fixing
  n = length(p)
  renum = 1:n
  shift = [2:n;1]
  w = flipdim(vertex(p),1)
  beta = 1 - flipdim(interiorangle(p),1)
  while any(abs.(beta[n]-[0;1]).<eps()) & (n > 2)
    renum = renum[shift]
    w = w[shift]
    beta = beta[shift]
    if renum[1]==1
      deg = abs.(beta) .< eps()
      w[deg] = []
      beta[deg] = []
      renum = 1:2
      n = 2
      break
    end
  end

  p = Polygon(flipdim(w,1),1.-flipdim(beta,1))

  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)

  zeta, c = param(w,beta,zeta0,qdat)
  preprev = -c/abs(c)./flipdim(zeta,1)
  prevangle = angle.(preprev)

  # compute geometric stuff
  z = vertex(p)
  zmid = 0.5*(z+circshift(z,-1))
  dz = circshift(z,-1)-z
  area = 0.5*imag(sum(conj.(zmid).*dz))
  if area > 0
    Zc = -0.5im/area*sum(dz.*(abs.(zmid).^2+abs.(dz).^2/12))
  else
    Zc = mean(z)
  end
  J = 0.25*imag(sum(conj.(zmid).*dz.*(abs.(zmid).^2+abs.(dz).^2/12)))

  # first two entries are for c₁ and c₀.
  betaflip = flipdim(beta,1)
  ccoeff = Complex128[abs(c),0.0]
  mom = [sum(betaflip.*preprev)]
  for k = 1:ncoeff
    push!(mom,sum(betaflip.*preprev.^(k+1)))
    coeffk = abs(c)*getcoefflist(k+1,1,mom);
    push!(ccoeff,coeffk)
  end

  # fix a
  zetainf = 5.0*exp.(im*collect(0:π/2:3π/2))
  sigmainf = -c/abs(c)./zetainf
  zinf = evaluate_exterior(sigmainf,w,beta,zeta,c,qdat)
  ccoeff[2] = mean([zinf[i]-sum(ccoeff.*zetainf[i].^(1:-1:-ncoeff))
                  for i = 1:length(zetainf)])

  # first entry is d₀
  dcoeff = [dot(ccoeff,ccoeff)]
  for k = 1:ncoeff+1
    push!(dcoeff,dot(ccoeff[1:end-k],ccoeff[k+1:end]))
  end

  ExteriorMap(n,p.vert,p.angle,qdat,tol,zeta,c,preprev,prevangle,
              ncoeff,ccoeff,dcoeff,mom,area,Zc,J)
end

function Base.show(io::IO, m::ExteriorMap)
    println(io, "Exterior map with $(m.N) vertices")
    print(io,   "   ")
    print(io,"vertices: ")
    for i = 1:length(m.z)
        print(io,"($(round(real(m.z[i]),4)),$(round(imag(m.z[i]),4))), ")
    end
    println(io)
    #for i = 1:length(m.z)
    #println(io, "       $(round(m.z[i],4))")
    #end
    print(io,"   ")
    print(io, "interior angles/π: ")
    for i = 1:length(m.angle)
        print(io, "$(round(m.angle[i],4)), ")
    end
    println(io)
    print(io,"   ")
    print(io, "prevertices on circle: ")
    for i = length(m.ζ):-1:1
        print(io,"($(round(real(m.ζ[i]),4)),$(round(imag(m.ζ[i]),4))), ")
    end
    println(io)
    print(io, "   ")
    print(io, "prevertex angles/π: ")
    for i = 1:length(m.prevangle)
        print(io, "$(round(m.prevangle[i]/π,4)), ")
    end
    println(io)
    print(io, "   ")
    print(io, "constant = $(round(m.constant,4)), ")
    print(io, "accuracy = $(m.accuracy), ")
    println(io)
    print(io, "   ")
    print(io, "number of multipole coefficients = $(m.ncoeff)")
    println(io)

end

"""
    m::ExteriorMap(ζ::T[;inside::Bool=false]) -> T

Evaluates the mapping `m` at the vector of points `ζ`, which are
assumed to lie inside the unit circle if `inside` is `true`, or
are assumed outside the unit circle if `inside` is `false` (the default).

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> ζ = [0.1,0.5-0.75im,-0.25-0.3im];

julia> m(ζ;inside=true)
3-element Array{Complex{Float64},1}:
   -6.9344-7.68965im
 0.0439774-1.11249im
   2.41181-0.044779im

julia> ζ = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];

julia> m(ζ)
   3-element Array{Complex{Float64},1}:
      0.81614+3.02956im
     -2.25237-2.08523im
    -0.333104+0.975837im
```
"""
function (m::ExteriorMap)(ζ::Vector{Complex128};inside::Bool=false)
  if inside
    return evaluate_exterior(ζ,flipdim(m.z,1),1.-flipdim(m.angle,1),
            m.ζ,m.constant,m.qdata)
  else
    b = -m.constant/abs(m.constant)
    ζ[ζ.==0] = eps();
    ζ[abs.(ζ).<1] = ζ[abs.(ζ).<1]./abs.(ζ[abs.(ζ).<1])

    σ = b./ζ
    return evaluate_exterior(σ,flipdim(m.z,1),1.-flipdim(m.angle,1),
            m.ζ,m.constant,m.qdata)
  end

end

(m::ExteriorMap)(ζ::Complex128;inside::Bool=false) =
            getindex(m([ζ];inside=inside),1)



#= get various data about the map =#

"""
    length(m::ConformalMap) -> Integer

Returns the number of control points/vertices of the map `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> length(m)
4
```
"""
Base.length(m::ConformalMap) = m.N


"""
    parameters(m::ExteriorMap) -> Tuple{Vector{Complex128},Complex128}

Returns a tuple of a vector of the prevertices and the complex factor of
the exterior polygon mapping `m`.

# Example

```jldoctest paramtest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> prev, C = parameters(m);

julia> prev
4-element Array{Complex{Float64},1}:
       1.0+0.0im
  0.376406-0.926455im
 -0.902383-0.430935im
 -0.186756+0.982406im
```
"""
parameters(m::ExteriorMap) = flipdim(m.ζ,1), m.constant

"""
    coefficients(m::ConformalMap) -> Tuple{Vector{Complex128},Vector{Complex128}}

Returns a tuple of vectors of the complex coefficients of the multipole
expansion of the mapping \$z(\\zeta)\$ described by `m` as well as the
coefficients of the square magnitude of the mapping \$|z(\\zeta)|^2\$.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> ccoeff, dcoeff = coefficients(m);

julia> ccoeff
14-element Array{Complex{Float64},1}:
       1.0198+0.0im
    -0.210364-0.0161983im
  -0.00655708+0.0398156im
     0.136922+0.0951343im
    -0.095035+0.0891769im
    0.0184341+0.0299586im
    0.0136513+2.78095e-5im
   -0.0159533-0.00264418im
  -0.00167426-0.00501161im
  -0.00578705-0.000221652im
  -0.00447511+0.00252069im
   0.00469089-0.00150588im
  0.000441767-0.00192516im
 -0.000381357-0.00174291im
```
"""
coefficients(m::ConformalMap) = m.ccoeff, m.dcoeff

"""
    moments(m::ExteriorMap) -> Vector{Complex128}

Return the moments of the prevertices for exterior polygon mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> mom = moments(m)
13-element Array{Complex{Float64},1}:
 -2.46691e-9+3.04899e-9im
  -0.0128596+0.0780855im
    0.805587+0.559726im
    -1.12125+1.04835im
    0.316471+0.633964im
    0.462871+0.225702im
    -1.56266+0.0306815im
   -0.106975-0.476173im
   -0.720332-0.0496159im
     -1.1805+0.0838739im
      1.1618-0.762023im
  -0.0612155-0.5728im
   -0.223423-0.726949im
```
"""
moments(m::ExteriorMap) = m.mom

"""
    area(m::ConformalMap) -> Float64

Returns the area of the shape described by the mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> area(m)
2.9
```

```jldoctest
julia> c = Complex128[1];

julia> m = PowerMap(c);

julia> area(m)
3.141592653589793
```

"""
area(m::ConformalMap) = m.area


"""
    centroid(m::ConformalMap) -> Complex128

Returns the complex centroid position of the shape described by the
mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> centroid(m)
-0.20919540229885059 - 0.04022988505747128im
```
"""
centroid(m::ConformalMap) = m.Zc

"""
    Jmoment(m::ConformalMap) -> Float64

Returns the second area moment of the shape described by the
mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> Jmoment(m)
1.5768333333333333
```
"""
Jmoment(m::ConformalMap) = m.J

#= solving for Schwarz-Christoffel parameters (prevertices, constant) =#

function param(w::Vector{Complex128},beta::Vector{Float64},
                 zeta0::Vector{Complex128},
                 qdat::Tuple{Array{Float64,2},Array{Float64,2}})
  # w clockwise
  # beta turning angles

  n = length(w)
  if n == 2
    zeta = Complex128[-1,1]
  else
    len = abs.(diff(circshift(w,1)))
    nmlen = abs.(len[3:n-1]/len[2])
    if isempty(zeta0)
      y0 = zeros(n-1)
    else
      zeta0 = zeta0/zeta0[n]
      θ = angle.(zeta0)
      θ[θ.<=0] = θ[θ.<=0] + 2π
      dt = diff([0;θ[1:n-1];2π])
      @. y0 = log(dt[1:n-1]/dt[2:n])
    end

    depfun! = Depfun(beta,nmlen,qdat)

    F0 = similar(y0)
    df = OnceDifferentiable(depfun!, y0, F0)
    sol = nlsolve(df,y0,show_trace = :false)

    zeta = zeros(Complex128,n)
    θ = zeros(n-1)
    y_to_zeta!(zeta,θ,sol.zero)

  end

  mid = zeta[1]*exp(0.5*im*angle(zeta[2]/zeta[1]))
  dequad = DQuad(beta,qdat)
  c = (w[2]-w[1])/(dequad([zeta[1]],[mid],[1],zeta)[1]-dequad([zeta[2]],[mid],[2],zeta)[1])

  return zeta, c

end

function evaluate_exterior(zeta::Vector{Complex128},w::Vector{Complex128},
                  beta::Vector{Float64},prev::Vector{Complex128},
                  c::Complex128,qdat::Tuple{Array{Float64,2},Array{Float64,2}})

  # this assumes zeta inside the unit circle

  if isempty(zeta)
    nothing
  end

  n = length(w)
  neval = length(zeta)
  tol = 10.0^(-size(qdat[1],1))

  # set up the integrator
  dequad = DQuad(beta,qdat)

  # initialize the mapped evaluation points
  wp = zeros(Complex128,neval)

  # find the closest prevertices to each evaluation point and their
  #  corresponding distances
  dz = abs.(hcat([zeta for i=1:n]...)-vcat([transpose(prev) for i=1:neval]...))
  (dist,ind) = findmin(dz,2)
  sing = floor.(Int,(ind[:]-1)/neval)+1

  # find any prevertices in the evaluation list and set them equal to
  #  the corresponding vertices. The origin is also a singular point
  vertex = (dist[:] .< tol)
  wp[vertex] = w[sing[vertex]]
  zerop = abs.(zeta) .< tol
  wp[zerop] = Inf
  vertex = vertex .| zerop

  # the starting (closest) singularities for each evaluation point
  prevs = prev[sing]

  # set the initial values of the non-vertices
  wp[.!vertex] = w[sing[.!vertex]]

  # distance to singularity at the origin
  abszeta = abs.(zeta)

  # unfinished cases
  unf = .!vertex

  # set the integration starting points
  zetaold = copy(prevs)
  zetanew = copy(zetaold)
  dist = ones(neval)
  while any(unf)
    # distance to integrate still
    dist[unf] = min.(1,2*abszeta[unf]./abs.(zeta[unf]-zetaold[unf]))

    # new integration end point
    zetanew[unf] = zetaold[unf] + dist[unf].*(zeta[unf]-zetaold[unf])

    # integrate
    wp[unf] = wp[unf] + c*dequad(zetaold[unf],zetanew[unf],sing[unf],prev)

    # set new starting integration points for those that can be integrated
    #  further
    unf = dist .< 1
    zetaold[unf] = zetanew[unf]

    # only the first step can have a singularity
    sing .= 0

  end

  return wp

end

function evalderiv_exterior(zeta::Vector{Complex128},beta::Vector{Float64},
                  prev::Vector{Complex128},c::Complex128)

    n = length(prev)
    neval = length(zeta)
    beta = [beta;-2]
    terms = [hcat([1 .- zeta/prev[i] for i = 1:n]...) zeta]
    dz = c*exp.(log.(terms)*beta)
    terms2 = [[-beta[i]/prev[i] for i = 1:n];beta[n+1]]
    ddz = dz.*((1.0./terms)*terms2)
    return dz, ddz

end


"""
    evalderiv(zeta::Vector{Complex128},m::ExteriorMap,inside::Bool) -> Tuple{Vector{Complex128},Vector{Complex128}}

Evaluates the first and second derivatives of the mapping `m` at the vector
of points `zeta`, which are assumed to lie inside the unit circle if
`inside` is `true`, or are assumed outside the unit circle if `inside` is
`false`. The first entry in the tuple returned is the first derivative,
the second entry is the second derivative.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> zeta = [0.1,0.5-0.75im,-0.25-0.3im];

julia> dz, ddz = evalderiv(zeta,m,true);

julia> dz
3-element Array{Complex{Float64},1}:
  67.2068+76.6284im
 -1.11666+0.544576im
  3.99129-5.30641im
```
"""
function evalderiv(zeta::Vector{Complex128},m::ExteriorMap,inside::Bool)

  if inside
    return evalderiv(zeta,1.-flipdim(m.angle,1),m.ζ,m.constant)
  else
    b = -m.constant/abs(m.constant)
    zeta[zeta.==0] = eps();
    zeta[abs.(zeta).<1] = zeta[abs.(zeta).<1]./abs.(zeta[abs.(zeta).<1])

    sigma = b./zeta
    dsigma = -sigma./zeta
    ddsigma = -2.0*dsigma./zeta
    dz, ddz = evalderiv(sigma,1.-flipdim(m.angle,1),m.ζ,m.constant)
    ddz = ddz.*dsigma.^2 + dz.*ddsigma
    dz .*= dsigma
    return dz, ddz
  end
end

"""
    evalderiv(zeta::Vector{Complex128},m::ExteriorMap) -> Tuple{Vector{Complex128},Vector{Complex128}}

Evaluates the first and second derivatives of the mapping `m` at the vector
of points `zeta`, which are assumed to lie outside the unit circle.
The first entry in the tuple returned is the first derivative,
the second entry is the second derivative.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];

julia> dz, ddz = evalderiv(zeta,m);

julia> dz
3-element Array{Complex{Float64},1}:
 1.03053+0.00444992im
 1.00696-0.0115011im
 1.30078-0.266625im
```
"""
evalderiv(zeta::Vector{Complex128},m::ExteriorMap) = evalderiv(zeta,m,false)

"""
    evalderiv(zeta::Complex128,m...) -> Tuple{Complex128,Complex128}

Evaluates the derivatives of `m` at a single point `zeta`.
"""
function evalderiv(zeta::Complex128,m::ExteriorMap,x...)
  dz, ddz = evalderiv([zeta],m::ExteriorMap,x...)
  return dz[1],ddz[1]
end



function getcoefflist(power,div,mom)

  # mom are the moments, mom[1] is M₁, etc

  if power%div!=0
    error("Indivisible power")
  end
  pow = Int(power/div)

  # Find the set of multi-indices I for which
  # sum(k(t)*t) = power/div. Each row of I corresponds
  # to a different multi-index in the set
  I = [1]
  for j = 2:pow
    I = madvance(I)
  end

  # Find the coefficient for 1/zeta^(pow-1) in the fhat expansion
  coeff = 0
  for j = 1:size(I,1)
    sumI = sum(I[j,:])
    fact = 1
    for l = 1:pow
      il = I[j,l]
      fact *= mom[l]^il/factorial(il)/l^il
    end
    coeff += fact*(-1)^sumI
  end
  return -coeff/(pow-1)

end


function madvance(P)
  # Given a set P of multi-indices
  # that satisfy a certain moment condition
  # |p|_1 = m (for all p in P), find the set Pplus1
  # that satisfies |p'|_1 = m+1 (for all p' in Pplus1)
  nP = size(P,1)
  m = size(P,2)
  Ppad = [ones(nP) P zeros(nP)]
  rows = 0
  Pplus1 = Int64[]
  for k = 1:nP
    # Loop through each index in the kth
    #  member of P, and find non-zero indices for shifting
    for j = m:-1:0
      if Ppad[k,j+1] == 1
        rows += 1
        Pplus1 = [Pplus1;zeros(Int64,1,m+2)]
        Pplus1[rows,1:m+2] = Ppad[k,1:m+2]
        Pplus1[rows,j+2] = Pplus1[rows,j+2]+1
        Pplus1[rows,j+1] = Pplus1[rows,j+1]-1
      end
    end
  end
  Pplus1 = Pplus1[:,2:end]
  Pplus1 = sortrows(Pplus1)
  dP = sum(abs.([transpose(Pplus1[1,:]);diff(Pplus1,1)]),2)
  return Pplus1[dP[:].!=0,:]

end

#######

struct DabsQuad{T,N,NQ}
  beta :: Vector{T}
  qdat :: Tuple{Array{T,2},Array{T,2}}
end

function DabsQuad(beta::Vector{T},tol::T) where T
  n = length(beta)
  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)
  DabsQuad{T,n,nqpts}(beta,qdat)
end

function DabsQuad(beta::Vector{T},
                  qdat::Tuple{Array{T,2},Array{T,2}}) where T
  n = length(beta)
  nqpts = size(qdat[1],1)
  DabsQuad{T,n,nqpts}(beta,qdat)
end

function (I::DabsQuad{T,N,NQ})(zeta1::Vector{Complex128},zeta2::Vector{Complex128},sing1::Vector{Int64},zeta::Vector{Complex128}) where {T,N,NQ}

   (qnode,qwght) = I.qdat
   argz = angle.(zeta)

   argz1 = angle.(zeta1)
   argz2 = angle.(zeta2)
   ang21 = angle.(zeta2./zeta1)

   bigargz = transpose(argz[:,ones(Int,NQ)])

   discont = (argz2-argz1).*ang21 .< 0
   argz2[discont] += 2π*sign.(ang21[discont])

   if isempty(sing1)
     sing1 = zeros(size(zeta1))
   end
   result = zeros(Float64,size(zeta1))

   nontriv = find(zeta1.!=zeta2)
    #tic()
   for k in nontriv
     z1k, z2k, arg1k, arg2k, sing1k =
          zeta1[k], zeta2[k], argz1[k], argz2[k], sing1[k]
     zetas = vcat(zeta[1:sing1k-1],zeta[sing1k+1:end])
     dist = min(1,2*minimum(abs.(zetas-z1k))/abs(z2k-z1k))
     argr = arg1k + dist*(arg2k-arg1k)
     ind = ((sing1k+N) % (N+1)) + 1
     nd = 0.5*((argr-arg1k)*qnode[:,ind] + argr + arg1k)
     wt = 0.5*abs(argr-arg1k)*qwght[:,ind]
     θ = (nd[:,ones(Int,N)]-bigargz.+2π).%(2π)
     θ[θ.>π] = 2π-θ[θ.>π]
     terms = 2sin.(0.5θ)
     if !any(terms==0.0)
        if sing1k > 0
            terms[:,sing1k] ./= abs.(nd-arg1k)
            wt .*= (0.5*abs.(argr-arg1k)).^I.beta[sing1k]
        end
        result[k] = At_mul_B(exp.(log.(terms)*I.beta),wt)
        while dist < 1.0
            argl = argr
            zetal = exp(im*argl)
            dist = min(1,2*minimum(abs.(zeta-zetal)/abs(zetal-z2k)))
            argr = argl + dist*(arg2k-argl)
            nd = 0.5*((argr-argl)*qnode[:,N+1] + argr + argl)
            wt = 0.5*abs(argr-argl)*qwght[:,N+1]
            #θ = hcat([(nd - argz[i] + 2π).%(2π) for i = 1:N]...)
            θ = (nd[:,ones(Int,N)]-bigargz.+2π).%(2π)
            θ[θ.>π] = 2π-θ[θ.>π]
            terms = 2sin.(0.5θ)
            result[k] += At_mul_B(exp.(log.(terms)*I.beta),wt)
        end
     end
   end
    #toc()
    #nothing
   return result
end

struct DQuad{T,N,NQ}
  beta :: Vector{T}
  qdat :: Tuple{Array{T,2},Array{T,2}}
end

function DQuad(beta::Vector{T},tol::T) where T
  n = length(beta)
  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)
  DQuad{T,n,nqpts}(beta,qdat)
end

function DQuad(beta::Vector{T},
                  qdat::Tuple{Array{T,2},Array{T,2}}) where T
  n = length(beta)
  nqpts = size(qdat[1],1)
  DQuad{T,n,nqpts}(beta,qdat)
end


function (I::DQuad{T,N,NQ})(zeta1::Vector{Complex128},zeta2::Vector{Complex128},
          sing1::Vector{Int64},zeta::Vector{Complex128}) where {T,N,NQ}

   (qnode,qwght) = I.qdat

   beta = [I.beta;-2]

   bigzeta = transpose(zeta[:,ones(Int,NQ)])

   if isempty(sing1)
     sing1 = zeros(Int,size(zeta1))
   end
   result = zeros(Complex128,size(zeta1))

   nontriv = find(zeta1.!=zeta2)
    #tic()
   for k in nontriv
     z1k, z2k, sing1k = zeta1[k], zeta2[k], sing1[k]
     zetas = vcat(zeta[1:sing1k-1],zeta[sing1k+1:end])
     dist = min(1,2*minimum(abs.(zetas-z1k))/abs(z2k-z1k))
     zetar = z1k + dist*(z2k-z1k)
     ind = sing1k + (N+1)*(sing1k==0)

     nd = 0.5*((zetar-z1k)*qnode[:,ind] + zetar + z1k)
     wt = 0.5*(zetar-z1k)*qwght[:,ind]
     terms = 1 .- nd[:,ones(Int,N)]./bigzeta
     if !any(terms==0.0)
       terms = hcat(terms,nd)
        if sing1k > 0
            terms[:,sing1k] ./= abs.(terms[:,sing1k])
            wt .*= (0.5*abs.(zetar-z1k)).^beta[sing1k]
        end
        result[k] = At_mul_B(exp.(log.(terms)*beta),wt)

        while dist < 1.0
            zetal = zetar
            dist = min(1,2*minimum(abs.(zeta-zetal)/abs(zetal-z2k)))
            zetar = zetal + dist*(z2k-zetal)
            nd = 0.5*((zetar-zetal)*qnode[:,N+1] + zetar + zetal)
            wt = 0.5*(zetar-zetal)*qwght[:,N+1]
            terms = 1 .- nd[:,ones(Int,N)]./bigzeta
            terms = hcat(terms,nd)
            result[k] += At_mul_B(exp.(log.(terms)*beta),wt)
        end
      end
    end
   return result
end


struct Depfun{T,N,NQ}
    zeta :: Vector{Complex{T}}
    beta :: Vector{T}
    nmlen :: Vector{T}
    qdat :: Tuple{Array{T,2},Array{T,2}}
    θ    :: Vector{T}
    mid  :: Vector{Complex{T}}
    ints :: Vector{T}
    dabsquad :: DabsQuad{T,N,NQ}
end

function Depfun(beta::Vector{T},nmlen::Vector{T},qdat:: Tuple{Array{T,2},Array{T,2}}) where T
    # should compute nmlen in here
    n = length(beta)
    nqpts = size(qdat[1],1)
    zeta = zeros(Complex128,n)
    θ = zeros(n-1)
    mid = zeros(Complex128,n-2)
    ints = zeros(Complex128,n-2)
    dabsquad = DabsQuad(beta,qdat)
    Depfun{T,n,nqpts}(zeta,beta,nmlen,qdat,θ,mid,ints,dabsquad)
end

function (R::Depfun{T,N,NQ})(F,y) where {T,N,NQ}


  y_to_zeta!(R.zeta,R.θ,y)

  @. R.mid = exp(im*0.5*(R.θ[1:N-2]+R.θ[2:N-1]))

  #tic()
  R.ints .= R.dabsquad(R.zeta[1:N-2],R.mid,collect(1:N-2),R.zeta)
  R.ints .+= R.dabsquad(R.zeta[2:N-1],R.mid,collect(2:N-1),R.zeta)

  #toc()

  if N > 3
    @. F[1:N-3] = abs(R.ints[2:N-2])/abs(R.ints[1]) - R.nmlen
  end

  res = -sum(R.beta./R.zeta)/R.ints[1]
  @. F[N-2] = real(res)
  @. F[N-1] = imag(res)
end


function y_to_zeta!(zeta::Vector{Complex{T}},
                    θ::Vector{T},y::Vector{T}) where T

  cs = cumsum(cumprod([1;exp.(-y)]))
  n = length(cs)
  @. θ = 2π*cs[1:n-1]/cs[n]
  zeta[n] = 1.0
  @. zeta[1:n-1] = exp(im*θ)

end

end
