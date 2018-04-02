using Documenter, SchwarzChristoffel

makedocs(
    format =:html,
    sitename = "SchwarzChristoffel.jl",
    pages = [
        "Home" => "index.md",
        "Usage" => "usage.md",
        "Reference" => [
                   "polygons.md",
                   "exterior.md"
                   ]
    ],
    assets = ["assets/custom.css"],
    strict = true,
    doctest = false
)

deploydocs(
    deps = nothing,
    repo = "github.com/jdeldre/SchwarzChristoffel.jl.git",
    target = "build",
    make = nothing,
    julia = "0.6"
)
