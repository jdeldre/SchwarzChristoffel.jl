module SchwarzChristoffel

using Reexport

include("Properties.jl")
include("Integration.jl")
include("Polygons.jl")
@reexport using .Polygons

include("Exterior.jl")
@reexport using .Exterior

export Polygons, Exterior

include("Plotting.jl")

end
