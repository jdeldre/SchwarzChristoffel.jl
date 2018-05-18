# SchwarzChristoffel

| Documentation | Build Status |
|:---:|:---:|
| [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://github.com/jdeldre/SchwarzChristoffel.jl) | [![Build Status](https://img.shields.io/travis/jdeldre/SchwarzChristoffel.jl/master.svg?label=linux)](https://travis-ci.org/jdeldre/SchwarzChristoffel.jl) [![Build status](https://img.shields.io/appveyor/ci/jdeldre/SchwarzChristoffel-jl/master.svg?label=windows)](https://ci.appveyor.com/project/jdeldre/schwarzchristoffel-jl/branch/master) [![codecov](https://codecov.io/gh/jdeldre/SchwarzChristoffel.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jdeldre/SchwarzChristoffel.jl) |

## About the package

The purpose of this package is to enable easy construction and evaluation of the mapping from the region inside or outside the unit circle to the exterior of a closed polygon.

The engine for constructing the mapping and its inverse is based on the work of Driscoll and Trefethen, [Schwarz-Christoffel Mapping](http://www.math.udel.edu/~driscoll/research/conformal.html), Cambridge University Press, 2002.

## Installation

This package requires Julia `0.6-` and above.
It is a registered package, so it should be installed with:
```julia
julia> Pkg.add("SchwarzChristoffel")
```
Since it is still under development, you should run
```julia
julia> Pkg.update()
```
to get the most recent version of the library and its dependencies.

Examples can be found in the [documentation](https://jdeldre.github.io/SchwarzChristoffel.jl).
