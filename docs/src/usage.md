# Basic usage

```@meta
DocTestSetup = quote
srand(1)
end
```

```@setup mapconstruct
using SchwarzChristoffel
```

First, we create a polygon shape by specifying its vertices. Note that the vertices must be provided in counter-clockwise order.

```@repl mapconstruct
x = [-1.0,0.2,1.0,-1.0]; y = [-1.0,-1.0,0.5,1.0];
p = Polygon(x,y)
```

Let's plot the polygon to make sure it matches what we wanted.
```@repl mapconstruct
plot(p)
savefig("polygon4.svg",format="svg"); nothing # hide
```

```@raw html
<object data="polygon4.svg" type="image/svg+xml"></object>
```

Now, we create the map from the unit circle to the polygon.

```@repl mapconstruct
m = ExteriorMap(p)
```

Let's visualize what we've constructed. Here, we will inspect the
mapping from the exterior of the unit circle to the exterior of the polygon.

```@repl mapconstruct
conformal_grid(m)
savefig("polygongrid.svg",format="svg"); nothing # hide
```

```@raw html
<object data="polygongrid.svg" type="image/svg+xml"></object>
```

We can now easily evaluate the map at any place we like. It could be evaluated
outside the unit circle:
```@repl mapconstruct
zeta = 1.2 + 0.1im
evaluate(zeta,m)
```

or it could be evaluated inside the unit circle:
```@repl mapconstruct
zeta = 0.5 + 0.1im
evaluate(zeta,m,true)
```

We can also evaluate the first and second derivative of the map at any place(s).
Let's evaluate at a range of points outside the circle.
```@repl mapconstruct
zeta = collect(1.1:0.1:2.0) + 0.1im
dz,ddz = evalderiv(zeta,m);
dz
```
```@setup mapconstruct2
using SchwarzChristoffel
```

Now let's try a more interesting shape. Here's a star-shaped body
```@repl mapconstruct2
n = 8; dθ = 2π/(2n)
θ = collect(0:dθ:2π-dθ)
w = (1+0.3cos.(n*θ)).*exp.(im*θ)
p = Polygon(w)
plot(p)
savefig("polygon8.svg",format="svg"); nothing # hide
```
```@raw html
<object data="polygon8.svg" type="image/svg+xml"></object>
```

Construct the map and plot it
```@repl mapconstruct2
m = ExteriorMap(p)
conformal_grid(m)
savefig("polygongrid8.svg",format="svg"); nothing # hide
```
```@raw html
<object data="polygongrid8.svg" type="image/svg+xml"></object>
```
