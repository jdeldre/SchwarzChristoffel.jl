# Basic usage

```@meta
DocTestSetup = quote
srand(1)
end
```

```@setup mapconstruct
using SchwarzChristoffel
```

```@example mapconstruct
p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
plot(p)
savefig("polygon4.svg",format="svg"); nothing # hide
```

```@raw html
<object data="polygon4.svg" width="400" type="image/svg+xml"></object>
```
