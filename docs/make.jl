using Documenter, SchwarzChristoffel

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

makedocs(
    #format =:html,
    sitename = "SchwarzChristoffel.jl",
    pages = [
        "Home" => "index.md",
        "Basic Usage" => "usage.md",
        "Reference" => [
                   "polygons.md",
                   "exterior.md"
                   ]
    ],
    #assets = ["assets/custom.css"],
    clean = true,
    doctest = true,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    )
)

deploydocs(
     repo = "github.com/jdeldre/SchwarzChristoffel.jl.git",
     target = "build",
     deps = nothing,
     make = nothing
     #versions = "v^"
)
