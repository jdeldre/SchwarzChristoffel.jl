using Documenter, SchwarzChristoffel

makedocs(
    format =:html,
    sitename = "SchwarzChristoffel.jl",
    pages = [
        "Home" => "index.md"
    ],
    assets = ["assets/custom.css"],
    strict = true
)

if "DOCUMENTER_KEY" in keys(ENV)
  deploydocs(
    repo = "github.com/jdeldre/SchwarzChristoffel.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    julia = "0.6"
  )
end
