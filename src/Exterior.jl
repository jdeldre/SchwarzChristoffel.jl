module Exterior

using NLsolve
#using DifferentialEquations
using OrdinaryDiffEq # using this module to reduce slow precompile of full DifferentialEquations package

using Compat
using Compat: reverse, repeat, findmin
using Compat.LinearAlgebra
using Compat.Statistics

using ..MapTypes
using ..Polygons

include("Integration.jl")
using .Integration

include("Reindex.jl")
using .Reindex

export PowerSeries,PowerSeriesDerivatives,PowerMap,ExteriorMap,
        KarmanTrefftzMap,JoukowskiMap,
        summary,parameters,coefficients,
        moments,area,centroid,Jmoment,addedmass

include("exterior/powerseries.jl")
include("exterior/polygon.jl")


struct PowerMap <: ConformalMap
    "Power series"
    ps::PowerSeries

    "Derivatives of power series"
    dps::PowerSeriesDerivatives

    "number of plotting control points"
    N::Int

    "control point coordinates in circle space"
    ζ::Vector{ComplexF64}

    "control point coordinates in body-fixed space"
    z::Vector{ComplexF64}

    "map Jacobian in body-fixed coordinates"
    dzdζ::Vector{ComplexF64}

    "Area enclosed by the mapped shape"
    area      :: Float64

    "Centroid of the mapped shape"
    Zc        :: ComplexF64

    "2nd area moment of the mapped shape"
    J         :: Float64

    "Added mass tensor"
    Ma        :: Array{Float64,2}

end

circle(N) = [exp(im*2π*(i-1)/N) for i in 1:N]

@doc raw"""
    PowerMap(c::Vector{ComplexF64}[;N = 200]) <: ConformalMap

Create a power series map from the exterior of the unit
circle to the exterior of a shape defined by the power series coefficients `c`.

The form of the mapping is

```math
z(\zeta) = c_{1}\zeta + c_{0} + \sum_{j=1}^{N_{c}} \frac{c_{-j}}{\zeta^{j}}
```

The entries in `c` correspond as follows: `c[1]` ``\rightarrow c_{1}``,
`c[2]` ``\rightarrow c_{0}``, `c[3]` ``\rightarrow c_{-1}``, etc.

The resulting map `m` can be evaluated at a single or a vector of points `ζ`
with `m(ζ)`.

# Example

```jldoctest
julia> c = ComplexF64[1,0,1/4];

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
function PowerMap(ccoeff::Vector{ComplexF64}; N::Int = 200)
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
  dcoeff = [ccoeff'ccoeff]
  for k = 1:ncoeff+1
      push!(dcoeff,ccoeff[1:end-k]'ccoeff[k+1:end])
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

(minv::InverseMap{PowerMap})(z::Vector{ComplexF64}) =
        evalinv_exterior(z,minv.m.ps,minv.m.dps)

(minv::InverseMap{PowerMap})(z::ComplexF64) = getindex(minv([z]),1)


(dm::DerivativeMap{PowerMap})(ζ) = dm.m.dps(ζ)


function shape_moments(ps::PowerSeries)

  ncoeff = length(ps.ccoeff)-2
  k = -1:ncoeff
  l = -1:ncoeff
  kml = repeat(transpose(k), length(ps.ccoeff))-repeat(l, 1, length(ps.ccoeff))

  c = Reflect(ShiftReindex(ps.ccoeff,-2))
  d = Reflect(OddReindex(ps.dcoeff))

  area = -π*sum(k.*abs.(c(-k)).^2)

  if area > 0
    Zc = -π/area*sum(k.*c(-k).*d(k))
  else
    Zc = ComplexF64(0)
  end
  J = real(-0.5π*c(-k)'*d(-kml)*(l.*c(-l)))

  return area, Zc, J
end

#=  Karman-Trefftz map from unit circle to airfoil
  This is actually a composition of maps, first from the unit circle in the ζ
  plane to a circle of radius a at center ϵCexp(iδ) in the ζ̃ plane. The mapping is
  (z-νC)/(z+νC) = ((ζ̃-C)/(ζ̃+C))^ν. This gives an interior angle at the
  trailing edge (ζ̃ = C) of angle (2-ν)π. Note that ζ̃ = ϵCexp(iδ) + aζ.
  The basic Joukowski airfoil is specified by ν = 2.

=#
struct KarmanTrefftzMap <: ConformalMap

  "K-T power ν, generally <∼ 2"
  nu::Real

  "Displacement of circle from origin, as fraction of C. Should be smaller than ~0.2."
  epsilon::Real

  "Angular orientation of circle, π/2 <= δ <= π"
  delta::Real

  "Scaling factor C. Airfoil chord is approximately 4C"
  C::Real

  "number of plotting control points"
  N::Int

  "radius of circle in ζ̃ plane"
  a::Real

  "center of circle in ζ̃ plane"
  mu::Number

  "Power series"
  ps::PowerSeries

  "control point coordinates in circle space"
  ζ::Vector{ComplexF64}

  "control point coordinates in body-fixed space"
  z::Vector{ComplexF64}

  "map Jacobian in body-fixed coordinates"
  dzdζ::Vector{ComplexF64}

  "Area enclosed by the mapped shape"
  area      :: Float64

  "Centroid of the mapped shape"
  Zc        :: ComplexF64

  "2nd area moment of the mapped shape"
  J         :: Float64

  "Added mass tensor"
  Ma        :: Array{Float64,2}

end

@doc raw"""
    KarmanTrefftzMap(ν,ϵ,δ,C[;N = 200]) <: ConformalMap

Create a map from the exterior of the unit
circle to the exterior of a Karman-Trefftz airfoil.

The form of the mapping is

```math
\frac{z-\nu C}{z+\nu C}  =
\left(\frac{\tilde{\zeta}-C}{\tilde{\zeta}+C}\right)^\nu
```
where $\tilde{\zeta}$ are the coordinates in an intermediate plane, in which
the circle is of radius ``a`` and centered at ``\epsilon C e^{i\delta}``:
```math
\tilde{\zeta} = \epsilon C e^{i\delta} + a \zeta
```
Note that ``a/C \geq 1`` and is determined by the choices for ``\epsilon`` and ``\delta``.

The trailing edge angle, ``(2-\nu)\pi`` is specified by ``\nu``. The thickness
is controlled by ``\epsilon C \cos\delta`` and the camber by ``\epsilon C \sin\delta``.
The airfoil chord length is approximately ``4C``. Generally, ``\epsilon`` should be
much smaller than 1 and ``\delta`` between ``\pi/2`` and ``\pi``.

The resulting map `m` can be evaluated at a single or a vector of points `ζ`
with `m(ζ)`.

# Example

```jldoctest
julia> ν = 1.9; ϵ = 0.1; δ = π; C = 0.25;

julia> m = KarmanTrefftzMap(ν,ϵ,δ,C)
Karman-Trefftz map

julia> ζ = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];

julia> m(ζ)
3-element Array{Complex{Float64},1}:
   0.268188+0.764722im
  -0.624265-0.502634im
 -0.0390996+0.126737im
```
""" KarmanTrefftzMap
function KarmanTrefftzMap(ν,ϵ,δ,C; N::Int = 200)

  β = atan(ϵ*sin(δ)/(1-ϵ*cos(δ)))
  a_C = cos(β) + sqrt(ϵ^2-sin(β)^2)
  a = a_C*C
  μ = ϵ*C*exp(im*δ)

  # Compute power series coefficients
  ncoeff = 4
  ccoeff = zeros(ComplexF64,ncoeff)
  ccoeff[1] = a
  ccoeff[2] = μ
  ccoeff[3] = (ν^2-1)*C/(3*a_C)
  ccoeff[4] = -ccoeff[3]*μ/a

  # Coefficients of |z(ζ)|²
  # dcoeff[1] = d₀, dcoeff[2] = d₋₁, etc.
  # Note that d₁ = conj(d₋₁) = conj(dcoeff[2]), d₂ = conj(d₋₂) = conj(dcoeff[3])
  dcoeff = [dot(ccoeff,ccoeff)]
  for k = 1:ncoeff+1
      push!(dcoeff,dot(ccoeff[1:end-k],ccoeff[k+1:end]))
  end
  ps = PowerSeries(ccoeff,dcoeff)

  ζ = circle(N)
  ζ̃ = similar(ζ)
  z = similar(ζ)
  dzdζ = similar(ζ)

  @. ζ̃ = μ + a*ζ
  @. z = ν*C*((ζ̃+C)^ν+(ζ̃-C)^ν)/((ζ̃+C)^ν-(ζ̃-C)^ν)
  @. dzdζ = 4ν^2*C^2*a*(ζ̃^2-C^2)^(ν-1)/((ζ̃+C)^ν-(ζ̃-C)^ν)^2

  area, Zc, J = shape_moments(ps)

  Ma = addedmass(ps,area)

  KarmanTrefftzMap(ν,ϵ,δ,C,N,a,μ,ps,ζ,z,dzdζ,area, Zc, J, Ma)
end

JoukowskiMap(y...;x...) = KarmanTrefftzMap(2,y...;x...)


function (m::KarmanTrefftzMap)(ζ::Number)
  ζ̃ = m.mu + m.a*ζ
  return m.nu*m.C*((ζ̃+m.C)^m.nu+(ζ̃-m.C)^m.nu)/((ζ̃+m.C)^m.nu-(ζ̃-m.C)^m.nu)
end

function (minv::InverseMap{KarmanTrefftzMap})(z::Number)
  ζ̃ = minv.m.C*((z+minv.m.nu*minv.m.C)^(1/minv.m.nu)+(z-minv.m.nu*minv.m.C)^(1/minv.m.nu))/
           ((z+minv.m.nu*minv.m.C)^(1/minv.m.nu)-(z-minv.m.nu*minv.m.C)^(1/minv.m.nu))
  return  (ζ̃ - minv.m.mu)/minv.m.a
end

function (dm::DerivativeMap{KarmanTrefftzMap})(ζ::Number)
  ζ̃ = dm.m.mu + dm.m.a*ζ
  return 4dm.m.nu^2*dm.m.C^2*dm.m.a*(ζ̃^2-dm.m.C^2)^(dm.m.nu-1)/
                      ((ζ̃+dm.m.C)^dm.m.nu-(ζ̃-dm.m.C)^dm.m.nu)^2
end

(m::KarmanTrefftzMap)(ζ::Vector{T}) where T<:Number = m.(ζ)

(minv::InverseMap{KarmanTrefftzMap})(z::Vector{T}) where T<:Number = minv.(z)


(dm::DerivativeMap{KarmanTrefftzMap})(ζ::Vector{T}) where T<:Number = dm.(ζ)


function Base.show(io::IO, m::KarmanTrefftzMap)
    if m.nu == 2
      println(io, "Joukowski map")
    else
      println(io, "Karman-Trefftz map")
    end
end

function Base.summary(m::KarmanTrefftzMap)
  println("Karman-Trefftz map of unit circle to exterior of airfoil")
  print("   ")
  print("circle center: ")
  println("($(round(m.epsilon*cos(m.delta),4)),$(round(m.epsilon*sin(m.delta),4)))")
  print("   ")
  print("chord length: ")
  println("$(4*m.C)")
  print("   ")
  print("trailing edge angle (deg): ")
  println("$(round((2-m.nu)*180,4))")
end


#=   Exterior map from circle to polygon  =#

struct ExteriorMap <: ConformalMap

  "Number of vertices on polygon"
  N :: Int

  "Coordinates of vertices on the polygon, defined ccw"
  z :: Vector{ComplexF64}

  "Interior angles in the polygon"
  angle  :: Vector{Float64}

  "Quadrature nodes and weights for the Gauss-Jacobi"
  qdata  :: Tuple{Array{Float64,2},Array{Float64,2}}

  "Accuracy of the quadrature"
  accuracy  :: Float64

  "Coordinates of the prevertices on the unit circle (interior)"
  ζ :: Vector{ComplexF64}

  "Constant factor of the mapping"
  constant  :: ComplexF64

  "Coordinates of the pre-prevertices on the unit circle (exterior)"
  preprev   :: Vector{ComplexF64}

  "Angular locations of the pre-prevertices on the unit circle (exterior)"
  prevangle :: Vector{Float64}

  "Number of multipole coefficients"
  ncoeff :: Int

  "Multipole coefficients for z(ζ) and |z(ζ)|²"
  ps    :: PowerSeries

  "Moments of prevertices"
  mom       :: Vector{ComplexF64}

  "Area enclosed by the mapped polygon"
  area      :: Float64

  "Centroid of the mapped polygon"
  Zc        :: ComplexF64

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

`ExteriorMap(p;ncoeff=200)` manually sets the number of coefficients of
negative powers of the multipole expansion of the mapping to `200`
(the default is 100).

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
function ExteriorMap(p::Polygon;tol::Float64 = 1e-8,ncoeff::Int = 100)

  n = length(p.vert)

  zeta0 = ComplexF64[]

  w = reverse(vertex(p), dims = 1)
  beta = 1 .- reverse(interiorangle(p), dims = 1)

  # do some fixing
  n = length(p)
  renum = 1:n
  shift = [2:n;1]
  w = reverse(vertex(p), dims = 1)
  beta = 1 .- reverse(interiorangle(p), dims = 1)
  while any(abs.(beta[n] .- [0;1]) .< eps()) & (n > 2)
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

  p = Polygon(reverse(w, dims = 1),1 .- reverse(beta, dims = 1))

  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)

  zeta, c = param(w,beta,zeta0,qdat)
  preprev = -c/abs(c)./reverse(zeta, dims = 1)
  prevangle = angle.(preprev)

  # first two entries are for c₁ and c₀.
  betaflip = reverse(beta, dims = 1)
  mom = [sum(betaflip.*preprev)]
  for k = 1:ncoeff
    push!(mom,sum(betaflip.*preprev.^(k+1)))
  end

  # use numerical quadrature to compute the first two coefficients
  # the first coefficient should be abs(c), but we don't set this
  # explicitly
  nθ = 1024
  dθ = 2π/nθ
  ζs = exp.(im*collect(0:nθ-1)*dθ)
  σs = -c/abs(c)./ζs
  dc = evaluate_exterior(σs,w,beta,zeta,c,qdat).*sin(dθ)
  dc ./= ζs
  ccoeff = ComplexF64[]
  for k = -1:0
    push!(ccoeff,sum(dc))
    dc .*= ζs
  end
  ccoeff ./= 2π

  # now, we use Gauss-Jacobi quadrature for the remaining coefficients.
  ζ = [zeta;zeta[1]]
  N = n+1
  argz1 = angle.(ζ[1:N-1])
  argz2 = angle.(ζ[2:N])
  ang21 = angle.(ζ[2:N]./ζ[1:N-1])
  discont = (argz2-argz1).*ang21 .< 0
  argz2[discont] += 2π*sign.(ang21[discont])
  mid = exp.(im*0.5*(argz1+argz2))
  qdat2 = qdata(beta,max(ceil(Int,ncoeff/2),nqpts))
  dabsquad = DabsQuad(beta,qdat2)
  for k = 1:ncoeff
     push!(ccoeff,c*(-c/abs(c))^k*sum(dabsquad(ζ[1:N-1],mid,collect(1:N-1),zeta,-k)+
                  dabsquad(ζ[2:N],mid,[collect(2:N-1);1],zeta,-k))/(2π*k))
  end

  # Now compute the coefficients of the square mapping
  # first entry is d₀
  dcoeff = [ccoeff'ccoeff]
  for k = 1:ncoeff+1
    push!(dcoeff,ccoeff[1:end-k]'ccoeff[k+1:end])
  end
  ps = PowerSeries(ccoeff,dcoeff)

  area, Zc, J = shape_moments(vertex(p))

  Ma = addedmass(ps,area)

  ExteriorMap(n,p.vert,p.angle,qdat,tol,zeta,c,preprev,prevangle,
              ncoeff,ps,mom,area,Zc,J,Ma)
end

function shape_moments(z::Vector{ComplexF64})

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
    print("($(round(real(m.z[i]), digits=4)),$(round(imag(m.z[i]), digits=4))), ")
  end
  println()
  print("   ")
  print("interior angles/π: ")
  for i = 1:length(m.angle)
    print("$(round(m.angle[i], digits=4)), ")
  end
  println()
  print("   ")
  print("prevertices on circle: ")
  for i = length(m.ζ):-1:1
    print("($(round(real(m.ζ[i]), digits=4)),$(round(imag(m.ζ[i]), digits=4))), ")
  end
  println()
  print("   ")
  print("prevertex angles/π: ")
  for i = 1:length(m.prevangle)
    print("$(round(m.prevangle[i]/π, digits=4)), ")
  end
  println()
  print("   ")
  print("constant = $(round(m.constant, digits=4)), ")
  print("accuracy = $(m.accuracy), ")
  println()
  print("   ")
  print("number of multipole coefficients = $(m.ncoeff)")
  println()

end

function (m::ExteriorMap)(ζ::Vector{ComplexF64};inside::Bool=false)
  if inside
    return evaluate_exterior(ζ,reverse(m.z, dims = 1),1 .- reverse(m.angle, dims = 1),
            m.ζ,m.constant,m.qdata)
  else
    b = -m.constant/abs(m.constant)
    ζ[ζ.==0] .= eps();
    ζ[abs.(ζ).<1] = ζ[abs.(ζ).<1]./abs.(ζ[abs.(ζ).<1])

    σ = b./ζ
    return evaluate_exterior(σ,reverse(m.z, dims = 1),1 .- reverse(m.angle, dims = 1),
            m.ζ,m.constant,m.qdata)
  end

end

(m::ExteriorMap)(ζ::ComplexF64;inside::Bool=false) = getindex(m([ζ];inside=inside),1)

"""
    InverseMap(m::ConformalMap)

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

function (minv::InverseMap{ExteriorMap})(z::Vector{ComplexF64};inside::Bool=false)
  if inside
    return evalinv_exterior(z,reverse(minv.m.z, dims = 1),1 .- reverse(minv.m.angle, dims = 1),
            minv.m.ζ,minv.m.constant,minv.m.qdata)
  else
    σ = evalinv_exterior(z,reverse(minv.m.z, dims = 1),1 .- reverse(minv.m.angle, dims = 1),
            minv.m.ζ,minv.m.constant,minv.m.qdata)
    b = -minv.m.constant/abs(minv.m.constant)

    return b./σ
  end

end

(minv::InverseMap{ExteriorMap})(z::ComplexF64;inside::Bool=false) =
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
```
"""
function DerivativeMap() end

function (dm::DerivativeMap{ExteriorMap})(ζ::Vector{ComplexF64};inside::Bool=false)
  if inside
    return evalderiv_exterior(ζ,1 .- reverse(dm.m.angle, dims = 1),dm.m.ζ,dm.m.constant)
  else
    b = -dm.m.constant/abs(dm.m.constant)
    ζ[ζ.==0] .= eps();
    ζ[abs.(ζ).<1] = ζ[abs.(ζ).<1]./abs.(ζ[abs.(ζ).<1])

    σ = b./ζ
    dσ = -σ./ζ
    ddσ = -2.0*dσ./ζ
    dz, ddz = evalderiv_exterior(σ,1 .- reverse(dm.m.angle, dims = 1),dm.m.ζ,dm.m.constant)
    ddz = ddz.*dσ.^2 + dz.*ddσ
    dz .*= dσ
    return dz, ddz
  end
end

function (dm::DerivativeMap{ExteriorMap})(ζ::ComplexF64;inside::Bool=false)
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
   number of multipole coefficients = 100
```

"""
Base.summary


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
    parameters(m::ExteriorMap) -> Tuple{Vector{ComplexF64},ComplexF64}

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
parameters(m::ExteriorMap) = reverse(m.ζ, dims = 1), m.constant

"""
    coefficients(m::ConformalMap) -> Tuple{Vector{ComplexF64},Vector{ComplexF64}}

Returns a tuple of vectors of the complex coefficients of the multipole
expansion of the mapping \$z(\\zeta)\$ described by `m` as well as the
coefficients of the square magnitude of the mapping \$|z(\\zeta)|^2\$.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> ccoeff, dcoeff = coefficients(m);
```
"""
coefficients(m::ConformalMap) = m.ps.ccoeff, m.ps.dcoeff

"""
    moments(m::ExteriorMap) -> Vector{ComplexF64}

Return the moments of the prevertices for exterior polygon mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> mom = moments(m);
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
julia> c = ComplexF64[1];

julia> m = PowerMap(c);

julia> area(m)
3.141592653589793
```

"""
area(m::ConformalMap) = m.area


"""
    centroid(m::ConformalMap) -> ComplexF64

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
  0.725129    0.0944902  -1.37387
  0.0944902   3.67634    -0.255119
 -1.37387    -0.255119    3.59231
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


end
