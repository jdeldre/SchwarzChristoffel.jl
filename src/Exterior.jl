module Exterior

using NLsolve

using ..Properties
using ..Polygons

include("Integration.jl")
using .Integration

include("Reindex.jl")
using .Reindex

export PowerSeries,PowerSeriesDerivatives,PowerMap,ExteriorMap,
        summary,parameters,coefficients,
        moments,area,centroid,Jmoment,addedmass

struct PowerSeries{N,T}
  ccoeff :: Vector{T}
  dcoeff :: Vector{T}
end

function PowerSeries(ccoeff::Vector{T},dcoeff::Vector{T}) where T
  ncoeff = length(ccoeff)-2
  PowerSeries{ncoeff,T}(ccoeff,dcoeff)
end

function Base.show(io::IO,ps::PowerSeries{N,T}) where {N,T}
  println(io, "multipole coefficients:")
  println(io, "  c₁ = $(ps.ccoeff[1]), ")
  println(io, "  c₀ = $(ps.ccoeff[2]), ")
  print(io,   "  c₋ᵢ = ")
  for i = 1:N
    print(io,"$(ps.ccoeff[2+i]), ")
  end
  println(io, "i = 1:$(N)")
end

function (ps::PowerSeries)(ζ::T) where T<:Number
  ζⁿ = ζ
  z = zero(ζ)
  for c in ps.ccoeff
    z += c*ζⁿ
    ζⁿ /= ζ
  end
  return z
end

(ps::PowerSeries)(ζ::Vector{T}) where T<:Number = ps.(ζ)

struct PowerSeriesDerivatives
  ps :: PowerSeries
end

function (dps::PowerSeriesDerivatives)(ζ::T) where T<:Number
  C = dps.ps.ccoeff
  dz = C[1]
  ζⁿ = 1/ζ^2
  ddz = Complex128(0)
  for n in 1:length(C)-2
    dz -= n*C[n+2]*ζⁿ
    ζⁿ /= ζ
    ddz += n*(n+1)*C[n+2]*ζⁿ
  end
  return dz, ddz
end

function (dps::PowerSeriesDerivatives)(ζs::Vector{T}) where T<:Number
  dz = zeros(ζs)
  ddz = zeros(ζs)
  for (i,ζ) in enumerate(ζs)
    dz[i], ddz[i] = dps(ζ)
  end
  return dz, ddz
end


struct PowerMap <: ConformalMap
    "Power series"
    ps::PowerSeries

    "Derivatives of power series"
    dps::PowerSeriesDerivatives

    "number of plotting control points"
    N::Int

    "control point coordinates in circle space"
    ζ::Vector{Complex128}

    "control point coordinates in body-fixed space"
    z::Vector{Complex128}

    "map Jacobian in body-fixed coordinates"
    dzdζ::Vector{Complex128}

    "Area enclosed by the mapped shape"
    area      :: Float64

    "Centroid of the mapped shape"
    Zc        :: Complex128

    "2nd area moment of the mapped shape"
    J         :: Float64

    "Added mass tensor"
    Ma        :: Array{Float64,2}

end

circle(N) = [exp(im*2π*(i-1)/N) for i in 1:N]

doc"""
    PowerMap(c::Vector{Complex12}[;N = 200]) <: ConformalMap

Create a power series map from the exterior of the unit
circle to the exterior of a shape defined by the power series coefficients `c`.

The form of the mapping is

```math
z(\zeta) = c_{1}\zeta + c_{0} + \sum_{j=1}^{N_{c}} \frac{c_{-j}}{\zeta^{j}}
```

The entries in `c` correspond as follows: `c[1]`$\rightarrow c_{1}$,
`c[2]`$\rightarrow c_{0}$, `c[3]`$\rightarrow c_{-1}$, etc.

The resulting map `m` can be evaluated at a single or a vector of points `ζ`
with `m(ζ)`.

# Example

```jldoctest
julia> c = Complex128[1,0,1/4];

julia> m = PowerMap(c)
Power series map

julia> ζ = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];

julia> m(ζ)
3-element Array{Complex{Float64},1}:
   1.025+2.925im
 -2.0625-1.9375im
     0.0+0.872727im
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


  # Coefficients of |z(ζ)|²
  # dcoeff[1] = d₀, dcoeff[2] = d₋₁, etc.
  # Note that d₁ = conj(d₋₁) = conj(dcoeff[2]), d₂ = conj(d₋₂) = conj(dcoeff[3])
  dcoeff = [dot(ccoeff,ccoeff)]
  for k = 1:ncoeff+1
      push!(dcoeff,dot(ccoeff[1:end-k],ccoeff[k+1:end]))
  end
  ps = PowerSeries(ccoeff,dcoeff)
  dps = PowerSeriesDerivatives(ps)

  ζ = circle(N)
  z = ps(ζ)
  dz, ddz = dps(ζ)


  area, Zc, J = shape_moments(ps)

  Ma = addedmass(ps,area)

  PowerMap(ps, dps, N, ζ, z, dz, area, Zc, J, Ma)
end

function Base.show(io::IO, m::PowerMap)
    println(io, "Power series map")
end

function Base.summary(m::PowerMap)
  m.ps
end

(m::PowerMap)(ζ) = m.ps(ζ)

(dm::DerivativeMap{PowerMap})(ζ) = dm.m.dps(ζ)


function shape_moments(ps::PowerSeries)

  ncoeff = length(ps.ccoeff)-2
  k = -1:ncoeff
  l = -1:ncoeff
  kml = k[:,ones(Int,length(ps.ccoeff))]'-l[:,ones(Int,length(ps.ccoeff))]

  c = Reflect(ShiftReindex(ps.ccoeff,-2))
  d = Reflect(OddReindex(ps.dcoeff))

  area = -π*sum(k.*abs.(c(-k)).^2)

  if area > 0
    Zc = -π/area*sum(k.*c(-k).*d(k))
  else
    Zc = Complex128(0)
  end
  J = Float64(-0.5π*c(-k)'*d(-kml)*(l.*c(-l)))

  return area, Zc, J
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
  ncoeff :: Int

  "Multipole coefficients for z(ζ) and |z(ζ)|²"
  ps    :: PowerSeries

  "Moments of prevertices"
  mom       :: Vector{Complex128}

  "Area enclosed by the mapped polygon"
  area      :: Float64

  "Centroid of the mapped polygon"
  Zc        :: Complex128

  "2nd area moment of the mapped polygon"
  J         :: Float64

  "Added mass tensor"
  Ma        :: Array{Float64,2}

end

"""
    ExteriorMap(p::Polygon[;tol::Float64][,ncoeff::Int]) <: ConformalMap

Create a Schwarz-Christoffel map from the interior or exterior of
the unit circle to the exterior of polygon `p`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p)
Schwarz-Christoffel map of unit circle to exterior of polygon with 4 vertices
```

`ExteriorMap(p;tol=1e-12)` manually sets the tolerance to `1e-12`
(the default is 1e-8).

`ExteriorMap(p;ncoeff=16)` manually sets the number of coefficients of
negative powers of the multipole expansion of the mapping to `16`
(the default is 12).

The resulting map `m` can be evaluated at a single or vector of points `ζ`
with `m(ζ[;inside::Bool])`. The points are assumed to lie outside the unit
circle, unless the optional argument `inside=true`, in which case they are
assumed to lie inside the circle.

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
  ps = PowerSeries(ccoeff,dcoeff)

  area, Zc, J = shape_moments(vertex(p))

  Ma = addedmass(ps,area)

  ExteriorMap(n,p.vert,p.angle,qdat,tol,zeta,c,preprev,prevangle,
              ncoeff,ps,mom,area,Zc,J,Ma)
end

function shape_moments(z::Vector{Complex128})

  zmid = 0.5*(z+circshift(z,-1))
  dz = circshift(z,-1)-z
  area = 0.5*imag(sum(conj.(zmid).*dz))
  if area > 0
    Zc = -0.5im/area*sum(dz.*(abs.(zmid).^2+abs.(dz).^2/12))
  else
    Zc = mean(z)
  end
  J = 0.25*imag(sum(conj.(zmid).*dz.*(abs.(zmid).^2+abs.(dz).^2/12)))

  return area, Zc, J
end

function Base.show(io::IO, m::ExteriorMap)
    println(io, "Schwarz-Christoffel map of unit circle to exterior of"*
                " polygon with $(m.N) vertices")
end

function Base.summary(m::ExteriorMap)

  println("Schwarz-Christoffel map of unit circle to exterior"*
          " of polygon with $(m.N) vertices")
  print("   ")
  print("vertices: ")
  for i = 1:length(m.z)
    print("($(round(real(m.z[i]),4)),$(round(imag(m.z[i]),4))), ")
  end
  println()
  print("   ")
  print("interior angles/π: ")
  for i = 1:length(m.angle)
    print("$(round(m.angle[i],4)), ")
  end
  println()
  print("   ")
  print("prevertices on circle: ")
  for i = length(m.ζ):-1:1
    print("($(round(real(m.ζ[i]),4)),$(round(imag(m.ζ[i]),4))), ")
  end
  println()
  print("   ")
  print("prevertex angles/π: ")
  for i = 1:length(m.prevangle)
    print("$(round(m.prevangle[i]/π,4)), ")
  end
  println()
  print("   ")
  print("constant = $(round(m.constant,4)), ")
  print("accuracy = $(m.accuracy), ")
  println()
  print("   ")
  print("number of multipole coefficients = $(m.ncoeff)")
  println()

end

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

(m::ExteriorMap)(ζ::Complex128;inside::Bool=false) = getindex(m([ζ];inside=inside),1)

doc"""
    InverseMap(m::ExteriorMap)

Constructs the inverse conformal map of the conformal map `m`.

This inverse conformal map can be evaluated at a single or vector of points.
Points should be outside the body. Whether the resulting point in the circle
plane is interpreted inside or outside the circle is determined by the optional
argument `inside`, which defaults to `false`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> m⁻¹ = InverseMap(m);

julia> ζ = [1.0+3.0im,-2.0-2.0im,0.1+1.1im];

julia> m⁻¹(m(ζ))
3-element Array{Complex{Float64},1}:
  1.0+3.0im
 -2.0-2.0im
  0.1+1.1im
```
""" InverseMap

function (minv::InverseMap{ExteriorMap})(z::Vector{Complex128};inside::Bool=false)
  if inside
    return evalinv_exterior(z,flipdim(minv.m.z,1),1.-flipdim(minv.m.angle,1),
            minv.m.ζ,minv.m.constant,minv.m.qdata)
  else
    σ = evalinv_exterior(z,flipdim(minv.m.z,1),1.-flipdim(minv.m.angle,1),
            minv.m.ζ,minv.m.constant,minv.m.qdata)
    b = -minv.m.constant/abs(minv.m.constant)

    return b./σ
  end

end

(minv::InverseMap{ExteriorMap})(z::Complex128;inside::Bool=false) =
                getindex(minv([z];inside=inside),1)


"""
    DerivativeMap(m::ConformalMap)

Constructs new conformal maps from the first and second derivatives of the
conformal map `m`.

These new conformal maps can be evaluated at a single or vector of points just as
 `m` is. The first entry in the tuple returned is the first derivative,
the second entry is the second derivative.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> dm = DerivativeMap(m);

julia> ζ = [0.1,0.5-0.75im,-0.25-0.3im];

julia> dz, ddz = dm(ζ;inside=true);

julia> dz
3-element Array{Complex{Float64},1}:
  67.2068+76.6284im
 -1.11666+0.544576im
  3.99129-5.30641im

julia> ζ = 1.0+3.0im;

julia> dz, ddz = dm(ζ);

julia> dz
1.0305280030434558 + 0.0044499240190600114im
```
"""function DerivativeMap() end

function (dm::DerivativeMap{ExteriorMap})(ζ::Vector{Complex128};inside::Bool=false)
  if inside
    return evalderiv_exterior(ζ,1.-flipdim(dm.m.angle,1),dm.m.ζ,dm.m.constant)
  else
    b = -dm.m.constant/abs(dm.m.constant)
    ζ[ζ.==0] = eps();
    ζ[abs.(ζ).<1] = ζ[abs.(ζ).<1]./abs.(ζ[abs.(ζ).<1])

    σ = b./ζ
    dσ = -σ./ζ
    ddσ = -2.0*dσ./ζ
    dz, ddz = evalderiv_exterior(σ,1.-flipdim(dm.m.angle,1),dm.m.ζ,dm.m.constant)
    ddz = ddz.*dσ.^2 + dz.*ddσ
    dz .*= dσ
    return dz, ddz
  end
end

function (dm::DerivativeMap{ExteriorMap})(ζ::Complex128;inside::Bool=false)
  dz, ddz = dm([ζ];inside=inside)
  return dz[1],ddz[1]
end


#= get various data about the map =#
"""
    summary(m::ConformalMap)

Returns a summary of data for a conformal map

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> summary(m)
Schwarz-Christoffel map of unit circle to exterior of polygon with 4 vertices
   vertices: (-1.0,-1.0), (0.2,-1.0), (1.0,0.5), (-1.0,1.0),
   interior angles/π: 0.5, 0.656, 0.422, 0.422,
   prevertices on circle: (1.0,0.0), (0.3764,-0.9265), (-0.9024,-0.4309), (-0.1868,0.9824),
   prevertex angles/π: -0.7291, -0.3519, 0.1291, 0.7111,
   constant = 0.6722 + 0.7669im, accuracy = 1.0e-8,
   number of multipole coefficients = 12
```

"""function Base.summary() end


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
coefficients(m::ConformalMap) = m.ps.ccoeff, m.ps.dcoeff

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

"""
    addedmass(m::ConformalMap) -> Array{Float64,2}

Returns the added mass matrix of the shape described by the
conformal mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> addedmass(m)
3×3 Array{Float64,2}:
  0.723706    0.0944715  -1.3739
  0.0944715   3.67642    -0.255121
 -1.3739     -0.255121    3.59239
```
"""
addedmass(m::ConformalMap) = m.Ma

function addedmass(ps::PowerSeries,area::Float64)
  M = zeros(3,3)
  c = Reflect(ShiftReindex(ps.ccoeff,-2))
  d = Reflect(OddReindex(ps.dcoeff))

  k = 1:length(ps.ccoeff)+1
  M[1,1] = 0.5sum(k.*abs.(d(-k)).^2)
  M[1,2] = M[2,1] = -imag(c(1)*d(-1))
  M[1,3] = M[3,1] =  real(c(1)*d(-1))
  M[2,3] = M[3,2] = -imag(c(1)*c(-1))
  M[2,2] = abs(c(1))^2-real(c(1)*c(-1))-area/(2π)
  M[3,3] = abs(c(1))^2+real(c(1)*c(-1))-area/(2π)
  return 2π*M

end

include("exterior/parameters.jl")

end
