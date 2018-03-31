# SchwarzChristoffel

| Documentation | Build Status |
|:---:|:---:|
| [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://github.com/jdeldre/SchwarzChristoffel.jl) | [![Build Status](https://img.shields.io/travis/jdeldre/SchwarzChristoffel.jl/master.svg?label=linux)](https://travis-ci.org/jdeldre/SchwarzChristoffel.jl) [![Build status](https://img.shields.io/appveyor/ci/jdeldre/SchwarzChristoffel-jl/master.svg?label=windows)](https://ci.appveyor.com/project/jdeldre/schwarzchristoffel-jl/branch/master) [![codecov](https://codecov.io/gh/jdeldre/SchwarzChristoffel.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jdeldre/SchwarzChristoffel.jl) |

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

Examples can be found in the [documentation](https://darwindarak.github.io/PotentialFlow.jl) and the [Jupyter notebooks](https://github.com/darwindarak/PotentialFlow.jl/tree/master/examples).
