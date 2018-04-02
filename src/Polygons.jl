module Polygons


import Base: length, show, isinf

export Polygon,vertex,interiorangle,isinpoly

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
  vert :: Vector{Complex128}
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
Polygon(w::Vector{Complex128}) = Polygon(w,interiorangle(w))

"""
    vertex(p::Polygon) -> -> Vector{Complex128}

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



function show(io::IO, p::Polygon)
    println(io, "Polygon with $(length(p.vert)) vertices at")
    print(io,   "             ")
    for i = 1:length(p.vert)
        print(io,"($(real(p.vert[i])),$(imag(p.vert[i]))) ")
    end
    println(io)
    println(io, "             interior angles/π = $(round.(p.angle, 3))")
end

end
