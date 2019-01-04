module Polygons

using Compat
using Compat.Statistics: mean
using Compat: reverse
import Base: length, show, isinf

export Polygon,vertex,interiorangle,isinpoly,naca4

"""
    Polygon(x::Vector{Float64}, y::Vector{Float64})

A polygon defined by its vertices, which must be provided in
counter-clockwise order.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at
             (-1.0,-1.0) (0.2,-1.0) (1.0,0.5) (-1.0,1.0)
             interior angles/π = [0.5, 0.656, 0.422, 0.422]
```
"""
struct Polygon
  vert :: Vector{ComplexF64}
  angle :: Vector{Float64}

  Polygon(vert,angle) = abs(vert[end]-vert[1])<eps() ? new(vert[1:end-1],angle) : new(vert,angle)
end

Polygon(x::T,y::T,angle) where T = Polygon(x+im*y,angle)


Polygon(x::T,y::T) where T = Polygon(x+im*y,interiorangle(x+im*y))

"""
    Polygon(w::Vector{Complex128})

Sets up a polygon with the coordinates of the vertices specified
with complex vector `w`. As usual, these must be supplied in
counter-clockwise order.

# Example

```jldoctest
julia> p = Polygon([-1.0-1.0im,0.2-1.0im,1.0+0.5im,-1.0+1.0im])
Polygon with 4 vertices at
             (-1.0,-1.0) (0.2,-1.0) (1.0,0.5) (-1.0,1.0)
             interior angles/π = [0.5, 0.656, 0.422, 0.422]
```
"""
Polygon(w::Vector{ComplexF64}) = Polygon(w,interiorangle(w))

"""
    vertex(p::Polygon) -> Vector{Complex128}

Returns the vector of vertices of the polygon `p`, in complex form.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> vertex(p)
4-element Array{Complex{Float64},1}:
 -1.0-1.0im
  0.2-1.0im
  1.0+0.5im
 -1.0+1.0im
```
"""
vertex(p::Polygon) = p.vert

"""
    isinf(p::Polygon) -> Bool

Returns `true` if any vertex in polygon `p` is at infinity.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> isinf(p)
false
```
"""
Base.isinf(p::Polygon) = any(isinf.(vertex(p)))

"""
    length(p::Polygon) -> Integer

Returns the number of vertices of the polygon `p`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> length(p)
4
```
"""
Base.length(p::Polygon) = length(vertex(p))

"""
    interiorangle(p::Polygon) -> Vector{Float64}

Returns the vector of interior angles (divided by \$\\pi\$) of the polygon `p`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> interiorangle(p)
4-element Array{Float64,1}:
 0.5
 0.655958
 0.422021
 0.422021
```
"""
interiorangle(p::Polygon) = length(p.angle) != 0 ? p.angle : interiorangle(p.vertex)

function interiorangle(w::Vector{ComplexF64})
  if length(w)==0
    return []
  end

  atinf = isinf.(w)
  mask = .~(atinf .| circshift(atinf,-1) .| circshift(atinf,1))

  dw = w - circshift(w,1)
  dwshift = circshift(dw,-1)
  beta = fill(NaN,length(w))
  beta[mask] = mod.(angle.( -dw[mask].*conj.(dwshift[mask]) )/π,2)

  mods = abs.(beta .+ 1) .< 1e-12
  beta[mods] = fill!(similar(beta[mods]), 1)

  return beta

end


function isinpoly(z::ComplexF64,w::Vector{ComplexF64},beta::Vector{Float64},tol)

  index = 0.0

  scale = mean(abs.(diff(circshift(w,-1))))
  if ~any(scale > eps())
        return
  end
  w = w/scale
  z = z/scale
  d = w .- z
  d[abs.(d) .< eps()] .= eps()
  ang = angle.(circshift(d,-1)./d)/π
  tangents = sign.(circshift(w,-1)-w)

  # Search for repeated points and fix these tangents
  for p = findall( tangents .== 0 )
    v = [w[p+1:end];w]
    tangents[p] = sign(v[findfirst(v.!=w[p])]-w[p])
  end

  # Points close to an edge
  onbdy = abs.(imag.(d./tangents)) .< 10*tol
  onvtx = abs.(d) .< tol
  onbdy = onbdy .& ( (abs.(ang) .> 0.9) .| onvtx .| circshift(onvtx,-1) )

  if ~any(onbdy)
    index = round(sum(ang)/2)
  else
    S = sum(ang[.~onbdy])
    b = beta[onvtx]
    augment = sum(onbdy) - sum(onvtx) - sum(b)

    index = round(augment*sign(S) + S)/2
  end
  return index==1

end

isinpoly(z,w,beta) = isinpoly(z,w,beta,eps())

"""
    isinpoly(z::Complex128,p::Polygon) -> Bool

Returns `true` or `false` depending on whether `z` is inside
or outside polygon `p`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> isinpoly(0.0+0.0im,p)
true

julia> isinpoly(1.0+2.0im,p)
false
```
"""
isinpoly(z,p::Polygon) = isinpoly(z,p::Polygon,eps())

"""
    isinpoly(z::Complex128,p::Polygon,tol::Float64) -> Bool

Returns `true` if `z` is inside or within distance `tol` of polygon `p`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> isinpoly(-1.01+0.0im,p)
false

julia> isinpoly(-1.01+0.0im,p,1e-2)
true
```
"""
isinpoly(z,p::Polygon,tol) = isinpoly(z,p.vert,p.angle,tol)

winding(z,x...) = float.(isinpoly(z,x...))

#=  some specific shape families =#

"""
    naca4(cam,pos,t[;np=20][,Zc=0.0+0.0im][,len=1.0]) -> Vector{Complex128}

Generates the vertices of a NACA 4-digit airfoil of chord length 1. The
relative camber is specified by `cam`, the position of
maximum camber (as fraction of chord) by `pos`, and the relative thickness
by `t`.

The optional parameter `np` specifies the number of points on the upper
or lower surface. The optional parameter `Zc` specifies the mean position of
the vertices (which is set to the origin by default). The optional parameter
`len` specifies the chord length.

# Example

```jldoctest
julia> w = naca4(0.0,0.0,0.12);

julia> p = Polygon(w);
```
"""
function naca4(cam::Number,pos::Number,t::Number;np=20,Zc=0.0+0.0im,len=1.0)

# Here, cam is the fractional camber, pos is the fractional chordwise position
# of max camber, and t is the fractional thickness.

npan = 2*np-2

# Trailing edge bunching
an = 1.5
anp = an+1
x = zeros(np)

θ = zero(x)
yc = zero(x)

for j = 1:np
    frac = Float64((j-1)/(np-1))
    x[j] = 1 - anp*frac*(1-frac)^an-(1-frac)^anp;
    if x[j] < pos
        yc[j] = cam/pos^2*(2*pos*x[j]-x[j]^2)
        if pos > 0
            θ[j] = atan(2*cam/pos*(1-x[j]/pos))
        end
    else
        yc[j] = cam/(1-pos)^2*((1-2*pos)+2*pos*x[j]-x[j]^2)
        if pos > 0
            θ[j] = atan(2*cam*pos/(1-pos)^2*(1-x[j]/pos))
        end
    end
end

xu = zero(x)
yu = xu
xl = xu
yl = yu

yt = t/0.20*(0.29690*sqrt.(x)-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4)

xu = x-yt.*sin.(θ)
yu = yc+yt.*cos.(θ)

xl = x+yt.*sin.(θ)
yl = yc-yt.*cos.(θ)

rt = 1.1019*t^2;
θ0 = 0
if pos > 0
    θ0 = atan(2*cam/pos)
end
# Center of leading edge radius
xrc = rt*cos(θ0)
yrc = rt*sin(θ0)
θle = collect(0:π/50:2π)
xlec = xrc .+ rt*cos.(θle)
ylec = yrc .+ rt*sin.(θle)

# Assemble data
coords = [xu yu xl yl x yc]
cole = [xlec ylec]

# Close the trailing edge
xpanold = [0.5*(xl[np]+xu[np]); reverse(xl[2:np-1], dims = 1); xu[1:np-1]]
ypanold = [0.5*(yl[np]+yu[np]); reverse(yl[2:np-1], dims = 1); yu[1:np-1]]

xpan = zeros(npan)
ypan = zeros(npan)
for ipan = 1:npan
    if ipan < npan
        xpan1 = xpanold[ipan]
        ypan1 = ypanold[ipan]
        xpan2 = xpanold[ipan+1]
        ypan2 = ypanold[ipan+1]
    else
        xpan1 = xpanold[npan]
        ypan1 = ypanold[npan]
        xpan2 = xpanold[1]
        ypan2 = ypanold[1]
    end
    xpan[ipan] = 0.5*(xpan1+xpan2)
    ypan[ipan] = 0.5*(ypan1+ypan2)
end
w = ComplexF64[1;reverse(xpan, dims = 1)+im*reverse(ypan,dims = 1)]*len
w .-= sum(w)/length(w)
return w .+ Zc

end


function show(io::IO, p::Polygon)
    println(io, "Polygon with $(length(p.vert)) vertices at")
    print(io,   "             ")
    for i = 1:length(p.vert)
        print(io,"($(real(p.vert[i])),$(imag(p.vert[i]))) ")
    end
    println(io)
    println(io, "             interior angles/π = $(round.(p.angle, digits=3))")
end

end
