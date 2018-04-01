using Documenter, SchwarzChristoffel

makedocs(
    format =:html,
    sitename = "SchwarzChristoffel.jl",
    pages = [
        "Home" => ["index.md",
                   "reference.md"]
    ],
    assets = ["assets/custom.css"],
    strict = true
)

deploydocs(
    deps = nothing,
    repo = "github.com/jdeldre/SchwarzChristoffel.jl.git",
    target = "build",
    make = nothing,
    julia = "0.6"
)
