module SchwarzChristoffel



include("Properties.jl")
include("Integration.jl")
include("Polygons.jl")
using .Polygons

include("Exterior.jl")
using .Exterior

export Polygons



end
