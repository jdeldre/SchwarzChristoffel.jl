using Documenter, SchwarzChristoffel

makedocs(
    format =:html,
    sitename = "SchwarzChristoffel.jl",
    pages = [
        "Home" => "index.md",
        "Basic Usage" => "usage.md",
        "Reference" => [
                   "polygons.md",
                   "exterior.md"
                   ]
    ],
    assets = ["assets/custom.css"],
    strict = true,
    doctest = true
)

deploydocs(
    deps = nothing,
    repo = "github.com/jdeldre/SchwarzChristoffel.jl.git",
    target = "build",
    make = nothing,
    julia = "0.6"
)
