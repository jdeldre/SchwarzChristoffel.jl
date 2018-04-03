# SchwarzChristoffel

*A tool to map polygons.*


## About the package

The purpose of this package is to enable easy construction and evaluation of the conformal mapping from the region inside or outside the unit circle to the exterior of a closed polygon.

A polygon could be a simple shape, of course, like a square, with only a few vertices:
```@setup mapnaca
using SchwarzChristoffel
p = Polygon([-0.5,0.5,0.5,-0.5],[-0.5,-0.5,0.5,0.5])
m = ExteriorMap(p)
conformal_grid(m)
savefig("square.svg",format="svg")
```
![](square.svg)

or it could be a more complicated shape, like a NACA 4412 airfoil:
```@setup mapnaca
using SchwarzChristoffel
w = naca4(0.04,0.4,0.12;len=1)
p = Polygon(w)
m = ExteriorMap(p)
conformal_grid(m)
savefig("naca4412.svg",format="svg")
```
![](naca4412.svg)

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
