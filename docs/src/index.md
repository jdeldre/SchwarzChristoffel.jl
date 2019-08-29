# SchwarzChristoffel

*A tool to map polygons.*


## About the package

The purpose of this package is to enable easy construction and evaluation of the conformal mapping from the region inside or outside the unit circle to the exterior of a closed polygon.

A polygon could be a simple shape, of course, like a square, with only a few vertices:
```@setup mapnaca
using SchwarzChristoffel
using Plots
pyplot()
clibrary(:colorbrewer)
default(grid = false)
p = Polygon([-0.5,0.5,0.5,-0.5],[-0.5,-0.5,0.5,0.5])
m = ExteriorMap(p)
plot(m)
savefig("square.svg")
```
![](square.svg)

or it could be a more complicated shape, like a NACA 4412 airfoil (assembled with
  line segments between a finite number of points on its shape):
```@setup mapnaca
using SchwarzChristoffel
w = naca4(0.04,0.4,0.12;len=1)
p = Polygon(w)
m = ExteriorMap(p)
plot(m)
savefig("naca4412.svg")
```
![](naca4412.svg)

The engine for constructing the mapping and its inverse is based on the work of Driscoll and Trefethen, [Schwarz-Christoffel Mapping](http://www.math.udel.edu/~driscoll/research/conformal.html), Cambridge University Press, 2002.

## Installation

This package requires Julia `1.0` and above.
It is a registered package, so it should be installed with `add SchwarzChristoffel`
at the package manager prompt.

Examples can be found in the [documentation](https://jdeldre.github.io/SchwarzChristoffel.jl).
