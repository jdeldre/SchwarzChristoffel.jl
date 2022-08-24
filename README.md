# SchwarzChristoffel

| Documentation | Build Status |
|:---:|:---:|
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jdeldre.github.io/SchwarzChristoffel.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jdeldre.github.io/SchwarzChristoffel.jl/dev) | [![Build Status](https://github.com/jdeldre/SchwarzChristoffel.jl/workflows/CI/badge.svg)](https://github.com/jdeldre/SchwarzChristoffel.jl/actions) [![Coverage](https://codecov.io/gh/jdeldre/SchwarzChristoffel.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jdeldre/SchwarzChristoffel.jl) |

## About the package

The purpose of this package is to enable easy construction and evaluation of the mapping from the region inside or outside the unit circle to the exterior of a closed polygon.

The engine for constructing the mapping and its inverse is based on the work of Driscoll and Trefethen, [Schwarz-Christoffel Mapping](http://www.math.udel.edu/~driscoll/research/conformal.html), Cambridge University Press, 2002.

## Installation

This package requires Julia `1.0` and above.
It is a registered package, so it should be installed with `add SchwarzChristoffel`
at the package manager prompt.

Examples can be found in the [documentation](https://jdeldre.github.io/SchwarzChristoffel.jl).
