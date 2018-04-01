module Polygons

using PyPlot
import PyCall
import Base: length, show, isinf

export Polygon,vertex,interiorangle,isinpoly,plot

struct Polygon
  vert :: Vector{Complex128}
  angle :: Vector{Float64}

  Polygon(vert,angle) = abs(vert[end]-vert[1])<eps() ? new(vert[1:end-1],angle) : new(vert,angle)
end

Polygon(x::T,y::T,angle) where T = Polygon(x+im*y,angle)

"""
    Polygons.Polygon(x::Vector{Float64}, y::Vector{Float64})

Sets up a polygon with the coordinates of the vertices specified
with vectors `x` and `y`.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]
```
"""
Polygon(x::T,y::T) where T = Polygon(x+im*y,interiorangle(x+im*y))

"""
    Polygons.Polygon(w::Vector{Complex128})

Sets up a polygon with the coordinates of the vertices specified
with complex vector `w`.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0-1.0im,0.2-1.0im,1.0+0.5im,-1.0+1.0im])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]
```
"""
Polygon(w::Vector{Complex128}) = Polygon(w,interiorangle(w))

"""
    Polygons.vertex(p::Polygons.Polygon)

Returns the vector of vertices of the polygon `p`, in complex form.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]

julia> Polygons.vertex(p)
4-element Array{Complex{Float64},1}:
 -1.0-1.0im
  0.2-1.0im
  1.0+0.5im
 -1.0+1.0im
```
"""
vertex(p::Polygon) = p.vert

"""
    isinf(p::Polygons.Polygon)

Returns `true` if any vertex in polygon `p` is at infinity.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]

julia> isinf(p)
false
```
"""
Base.isinf(p::Polygon) = any(isinf.(vertex(p)))

"""
    length(p::Polygons.Polygon)

Returns the number of vertices of the polygon `p`.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]

julia> length(p)
4
```
"""
Base.length(p::Polygon) = length(vertex(p))

"""
    Polygons.interiorangle(p::Polygons.Polygon)

Returns the vector of interior angles of the polygon `p`.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]

julia> Polygons.interiorangle(p)
4-element Array{Float64,1}:
 0.5
 0.655958
 0.422021
 0.422021
```
"""
interiorangle(p::Polygon) = length(p.angle) != 0 ? p.angle : interiorangle(p.vertex)

function interiorangle(w::Vector{Complex128})
  if length(w)==0
    return []
  end

  atinf = isinf.(w)
  mask = .~(atinf .| circshift(atinf,-1) .| circshift(atinf,1))

  dw = w - circshift(w,1)
  dwshift = circshift(dw,-1)
  beta = fill(NaN,length(w))
  beta[mask] = mod.(angle.( -dw[mask].*conj.(dwshift[mask]) )/π,2)

  mods = abs.(beta+1) .< 1e-12
  beta[mods] = ones(beta[mods])

  return beta

end


function isinpoly(z::Complex128,w::Vector{Complex128},beta::Vector{Float64},tol)

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
  for p = find( tangents .== 0 )
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
    Polygons.isinpoly(z::Complex128,p::Polygons.Polygon)

Returns `true` or `false` depending on whether `z` is inside
or outside polygon `p`.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]

julia> Polygons.isinpoly(0.0+0.0im,p)
true

julia> Polygons.isinpoly(1.0+2.0im,p)
false
```

    Polygons.isinpoly(z::Complex128,p::Polygons.Polygon,tol::Float64)

Returns `true` if `z` is inside or within distance `tol` polygon `p`.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]

julia> Polygons.isinpoly(-1.01+0.0im,p)
false

julia> Polygons.isinpoly(-1.01+0.0im,p,1e-2)
true
```
"""
isinpoly(z,p::Polygon) = isinpoly(z,p::Polygon,eps())

isinpoly(z,p::Polygon,tol) = isinpoly(z,p.vert,p.angle,tol)

winding(z,x...) = float.(isinpoly(z,x...))





const mygreen = [151/255,180/255,118/255];
const mygreen2 = [113/255,161/255,103/255];
const myblue = [74/255,144/255,226/255];

PlotColor = Union{Vector{Float64},String,Symbol}

struct PlotStyle

  bodycolor :: PlotColor
  linewidth :: Int64
  pointsize :: Int64

end

function PlotStyle()
  PyCall.PyDict(matplotlib["rcParams"])["font.family"] = "STIXGeneral";
  PyCall.PyDict(matplotlib["rcParams"])["mathtext.fontset"] = "stix";
  PyCall.PyDict(matplotlib["rcParams"])["font.size"] = 10;

  bodycolor = mygreen
  linewidth = 1
  pointsize = 3

  PlotStyle(bodycolor,linewidth,pointsize)

end

"""
    Polygons.plot(p::Polygons.Polygon)

Plots the polygon `p`.

# Example

```jldoctest
julia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
Polygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]
             interior angles/π = [0.5, 0.656, 0.422, 0.422]

julia> Polygons.plot(p);
```
"""
function plot(p::Polygon)
  ps = PlotStyle()
  pf = PyPlot.fill(real(p.vert),imag(p.vert),facecolor=ps.bodycolor);
  PyPlot.plot(real(p.vert),imag(p.vert),color=ps.bodycolor);
  PyPlot.axis("scaled")
end

function show(io::IO, p::Polygon)
    println(io, "Polygon with $(length(p.vert)) vertices at $(p.vert) ")
    println(io, "             interior angles/π = $(round.(p.angle, 3))")
end

end
