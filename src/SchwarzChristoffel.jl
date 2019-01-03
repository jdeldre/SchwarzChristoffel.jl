module SchwarzChristoffel

using Reexport

include("MapTypes.jl")
@reexport using .MapTypes

include("Polygons.jl")
@reexport using .Polygons

include("Exterior.jl")
@reexport using .Exterior

export Polygons, Exterior

# plotting stuff
include("plot_recipes.jl")

end
