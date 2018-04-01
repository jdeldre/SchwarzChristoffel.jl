# SchwarzChristoffel

## About the package

The purpose of this package is to enable easy construction and evaluation of the mapping from the region inside or outside the unit circle to the exterior of a closed polygon.

## Installation

This package requires Julia `0.6-` and above.
It is not a registered package, so it should be installed with:
```julia
julia> Pkg.clone("git@github.com:jdeldre/SchwarzChristoffel.jl.git")
```
Since it is still under heavy development, you should run
```julia
julia> Pkg.test("SchwarzChristoffel") # might take some time
```
to make sure things are working as intended and
```julia
julia> Pkg.update()
```
to get the most recent version of the library and its dependencies.

Examples can be found in the [documentation](https://jdeldre.github.io/SchwarzChristoffel.jl) and the [Jupyter notebooks](https://github.com/jdeldre/SchwarzChristoffel.jl/tree/master/examples).
