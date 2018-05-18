var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SchwarzChristoffel-1",
    "page": "Home",
    "title": "SchwarzChristoffel",
    "category": "section",
    "text": "A tool to map polygons."
},

{
    "location": "index.html#About-the-package-1",
    "page": "Home",
    "title": "About the package",
    "category": "section",
    "text": "The purpose of this package is to enable easy construction and evaluation of the conformal mapping from the region inside or outside the unit circle to the exterior of a closed polygon.A polygon could be a simple shape, of course, like a square, with only a few vertices:using SchwarzChristoffel\nusing Plots\npyplot()\nclibrary(:colorbrewer)\ndefault(grid = false)\np = Polygon([-0.5,0.5,0.5,-0.5],[-0.5,-0.5,0.5,0.5])\nm = ExteriorMap(p)\nplot(m)\nsavefig(\"square.svg\")(Image: )or it could be a more complicated shape, like a NACA 4412 airfoil:using SchwarzChristoffel\nw = naca4(0.04,0.4,0.12;len=1)\np = Polygon(w)\nm = ExteriorMap(p)\nplot(m)\nsavefig(\"naca4412.svg\")(Image: )"
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia 0.6- and above. It is not a registered package, so it should be installed with:julia> Pkg.clone(\"git@github.com:jdeldre/SchwarzChristoffel.jl.git\")Since it is still under heavy development, you should runjulia> Pkg.test(\"SchwarzChristoffel\") # might take some timeto make sure things are working as intended andjulia> Pkg.update()to get the most recent version of the library and its dependencies.Examples can be found in the documentation and the Jupyter notebooks."
},

{
    "location": "usage.html#",
    "page": "Basic Usage",
    "title": "Basic Usage",
    "category": "page",
    "text": ""
},

{
    "location": "usage.html#Basic-usage-1",
    "page": "Basic Usage",
    "title": "Basic usage",
    "category": "section",
    "text": "DocTestSetup = quote\nsrand(1)\nendusing SchwarzChristoffel\nusing Plots\npyplot()\nclibrary(:colorbrewer)\ndefault(grid = false)First, we create a polygon shape by specifying its vertices. Note that the vertices must be provided in counter-clockwise order.x = [-1.0,0.2,1.0,-1.0]; y = [-1.0,-1.0,0.5,1.0];\np = Polygon(x,y)Let\'s plot the polygon to make sure it matches what we wanted.plot(p)\nsavefig(\"polygon4.svg\")(Image: )Now, we create the map from the unit circle to the polygon.m = ExteriorMap(p)Let\'s visualize what we\'ve constructed. Here, we will inspect the mapping from the exterior of the unit circle to the exterior of the polygon.plot(m)\nsavefig(\"polygongrid.svg\"); nothing # hide(Image: )We can now easily evaluate the map at any place we like. It could be evaluated outside the unit circle:ζ = 1.2 + 0.1im\nm(ζ)or it could be evaluated inside the unit circle:ζ = 0.5 + 0.1im\nm(ζ;inside=true)We can also evaluate the first and second derivative of the map at any place(s). Let\'s evaluate at a range of points outside the circle.dm = DerivativeMap(m)\nζ = collect(1.1:0.1:2.0) + 0.1im\ndz,ddz = dm(ζ);\ndzusing SchwarzChristoffel\nusing Plots\npyplot()\nclibrary(:colorbrewer)\ndefault(grid = false)Now let\'s try a more interesting shape. Here\'s a star-shaped bodyn = 8; dθ = 2π/(2n)\nθ = collect(0:dθ:2π-dθ)\nw = (1+0.3cos.(n*θ)).*exp.(im*θ)\np = Polygon(w)\nplot(p)\nsavefig(\"polygon8.svg\"); nothing # hide(Image: )Construct the map and plot itm = ExteriorMap(p)\nplot(m)\nsavefig(\"polygongrid8.svg\"); nothing # hide(Image: )"
},

{
    "location": "polygons.html#",
    "page": "Polygons",
    "title": "Polygons",
    "category": "page",
    "text": ""
},

{
    "location": "polygons.html#Polygons-1",
    "page": "Polygons",
    "title": "Polygons",
    "category": "section",
    "text": "DocTestSetup = quote\nusing SchwarzChristoffel\nsrand(1)\nend"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.Polygon",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.Polygon",
    "category": "type",
    "text": "Polygon(x::Vector{Float64}, y::Vector{Float64})\n\nA polygon defined by its vertices, which must be provided in counter-clockwise order.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at\n             (-1.0,-1.0) (0.2,-1.0) (1.0,0.5) (-1.0,1.0)\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\n\n\n"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.vertex",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.vertex",
    "category": "function",
    "text": "vertex(p::Polygon) -> Vector{Complex128}\n\nReturns the vector of vertices of the polygon p, in complex form.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> vertex(p)\n4-element Array{Complex{Float64},1}:\n -1.0-1.0im\n  0.2-1.0im\n  1.0+0.5im\n -1.0+1.0im\n\n\n\n"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.interiorangle",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.interiorangle",
    "category": "function",
    "text": "interiorangle(p::Polygon) -> Vector{Float64}\n\nReturns the vector of interior angles (divided by pi) of the polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> interiorangle(p)\n4-element Array{Float64,1}:\n 0.5\n 0.655958\n 0.422021\n 0.422021\n\n\n\n"
},

{
    "location": "polygons.html#Base.length-Tuple{SchwarzChristoffel.Polygons.Polygon}",
    "page": "Polygons",
    "title": "Base.length",
    "category": "method",
    "text": "length(p::Polygon) -> Integer\n\nReturns the number of vertices of the polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> length(p)\n4\n\n\n\n"
},

{
    "location": "polygons.html#Base.isinf-Tuple{SchwarzChristoffel.Polygons.Polygon}",
    "page": "Polygons",
    "title": "Base.isinf",
    "category": "method",
    "text": "isinf(p::Polygon) -> Bool\n\nReturns true if any vertex in polygon p is at infinity.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> isinf(p)\nfalse\n\n\n\n"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.isinpoly",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.isinpoly",
    "category": "function",
    "text": "isinpoly(z::Complex128,p::Polygon) -> Bool\n\nReturns true or false depending on whether z is inside or outside polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> isinpoly(0.0+0.0im,p)\ntrue\n\njulia> isinpoly(1.0+2.0im,p)\nfalse\n\n\n\nisinpoly(z::Complex128,p::Polygon,tol::Float64) -> Bool\n\nReturns true if z is inside or within distance tol of polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> isinpoly(-1.01+0.0im,p)\nfalse\n\njulia> isinpoly(-1.01+0.0im,p,1e-2)\ntrue\n\n\n\n"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.naca4",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.naca4",
    "category": "function",
    "text": "naca4(cam,pos,t[;np=20][,Zc=0.0+0.0im][,len=1.0]) -> Vector{Complex128}\n\nGenerates the vertices of a NACA 4-digit airfoil of chord length 1. The relative camber is specified by cam, the position of maximum camber (as fraction of chord) by pos, and the relative thickness by t.\n\nThe optional parameter np specifies the number of points on the upper or lower surface. The optional parameter Zc specifies the mean position of the vertices (which is set to the origin by default). The optional parameter len specifies the chord length.\n\nExample\n\njulia> w = naca4(0.0,0.0,0.12);\n\njulia> p = Polygon(w);\n\n\n\n"
},

{
    "location": "polygons.html#Methods-1",
    "page": "Polygons",
    "title": "Methods",
    "category": "section",
    "text": "Polygons.Polygon\nPolygons.vertex\nPolygons.interiorangle\nBase.length(::Polygons.Polygon)\nBase.isinf(::Polygons.Polygon)\nPolygons.isinpoly\nPolygons.naca4"
},

{
    "location": "polygons.html#Index-1",
    "page": "Polygons",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"polygons.md\"]"
},

{
    "location": "exterior.html#",
    "page": "Exterior map",
    "title": "Exterior map",
    "category": "page",
    "text": ""
},

{
    "location": "exterior.html#Exterior-map-1",
    "page": "Exterior map",
    "title": "Exterior map",
    "category": "section",
    "text": "DocTestSetup = quote\nusing SchwarzChristoffel\nsrand(1)\nend"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.ExteriorMap-Tuple{SchwarzChristoffel.Polygons.Polygon}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.ExteriorMap",
    "category": "method",
    "text": "ExteriorMap(p::Polygon[;tol::Float64][,ncoeff::Int]) <: ConformalMap\n\nCreate a Schwarz-Christoffel map from the interior or exterior of the unit circle to the exterior of polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p)\nSchwarz-Christoffel map of unit circle to exterior of polygon with 4 vertices\n\nExteriorMap(p;tol=1e-12) manually sets the tolerance to 1e-12 (the default is 1e-8).\n\nExteriorMap(p;ncoeff=200) manually sets the number of coefficients of negative powers of the multipole expansion of the mapping to 200 (the default is 100).\n\nThe resulting map m can be evaluated at a single or vector of points ζ with m(ζ[;inside::Bool]). The points are assumed to lie outside the unit circle, unless the optional argument inside=true, in which case they are assumed to lie inside the circle.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> ζ = [0.1,0.5-0.75im,-0.25-0.3im];\n\njulia> m(ζ;inside=true)\n3-element Array{Complex{Float64},1}:\n   -6.9344-7.68965im\n 0.0439774-1.11249im\n   2.41181-0.044779im\n\njulia> ζ = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];\n\njulia> m(ζ)\n3-element Array{Complex{Float64},1}:\n   0.81614+3.02956im\n  -2.25237-2.08523im\n -0.333104+0.975837im\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.PowerMap-Tuple{Array{Complex{Float64},1}}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.PowerMap",
    "category": "method",
    "text": "PowerMap(c::Vector{Complex12}[;N = 200]) <: ConformalMap\n\nCreate a power series map from the exterior of the unit circle to the exterior of a shape defined by the power series coefficients c.\n\nThe form of the mapping is\n\nz(zeta) = c_1zeta + c_0 + sum_j=1^N_c fracc_-jzeta^j\n\nThe entries in c correspond as follows: c[1]rightarrow c_1, c[2]rightarrow c_0, c[3]rightarrow c_-1, etc.\n\nThe resulting map m can be evaluated at a single or a vector of points ζ with m(ζ).\n\nExample\n\njulia> c = Complex128[1,0,1/4];\n\njulia> m = PowerMap(c)\nPower series map\n\njulia> ζ = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];\n\njulia> m(ζ)\n3-element Array{Complex{Float64},1}:\n   1.025+2.925im\n -2.0625-1.9375im\n     0.0+0.872727im\n\n\n\n"
},

{
    "location": "exterior.html#Base.summary-Tuple{}",
    "page": "Exterior map",
    "title": "Base.summary",
    "category": "method",
    "text": "summary(m::ConformalMap)\n\nReturns a summary of data for a conformal map\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> summary(m)\nSchwarz-Christoffel map of unit circle to exterior of polygon with 4 vertices\n   vertices: (-1.0,-1.0), (0.2,-1.0), (1.0,0.5), (-1.0,1.0),\n   interior angles/π: 0.5, 0.656, 0.422, 0.422,\n   prevertices on circle: (1.0,0.0), (0.3764,-0.9265), (-0.9024,-0.4309), (-0.1868,0.9824),\n   prevertex angles/π: -0.7291, -0.3519, 0.1291, 0.7111,\n   constant = 0.6722 + 0.7669im, accuracy = 1.0e-8,\n   number of multipole coefficients = 100\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.Jmoment-Tuple{SchwarzChristoffel.Properties.ConformalMap}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.Jmoment",
    "category": "method",
    "text": "Jmoment(m::ConformalMap) -> Float64\n\nReturns the second area moment of the shape described by the mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> Jmoment(m)\n1.5768333333333333\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.addedmass-Tuple{SchwarzChristoffel.Properties.ConformalMap}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.addedmass",
    "category": "method",
    "text": "addedmass(m::ConformalMap) -> Array{Float64,2}\n\nReturns the added mass matrix of the shape described by the conformal mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> addedmass(m)\n3×3 Array{Float64,2}:\n  0.725129    0.0944902  -1.37387\n  0.0944902   3.67634    -0.255119\n -1.37387    -0.255119    3.59231\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.area-Tuple{SchwarzChristoffel.Properties.ConformalMap}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.area",
    "category": "method",
    "text": "area(m::ConformalMap) -> Float64\n\nReturns the area of the shape described by the mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> area(m)\n2.9\n\njulia> c = Complex128[1];\n\njulia> m = PowerMap(c);\n\njulia> area(m)\n3.141592653589793\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.centroid-Tuple{SchwarzChristoffel.Properties.ConformalMap}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.centroid",
    "category": "method",
    "text": "centroid(m::ConformalMap) -> Complex128\n\nReturns the complex centroid position of the shape described by the mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> centroid(m)\n-0.20919540229885059 - 0.04022988505747128im\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.coefficients-Tuple{SchwarzChristoffel.Properties.ConformalMap}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.coefficients",
    "category": "method",
    "text": "coefficients(m::ConformalMap) -> Tuple{Vector{Complex128},Vector{Complex128}}\n\nReturns a tuple of vectors of the complex coefficients of the multipole expansion of the mapping z(zeta) described by m as well as the coefficients of the square magnitude of the mapping z(zeta)^2.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> ccoeff, dcoeff = coefficients(m);\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.moments-Tuple{SchwarzChristoffel.Exterior.ExteriorMap}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.moments",
    "category": "method",
    "text": "moments(m::ExteriorMap) -> Vector{Complex128}\n\nReturn the moments of the prevertices for exterior polygon mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> mom = moments(m);\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.parameters-Tuple{SchwarzChristoffel.Exterior.ExteriorMap}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.parameters",
    "category": "method",
    "text": "parameters(m::ExteriorMap) -> Tuple{Vector{Complex128},Complex128}\n\nReturns a tuple of a vector of the prevertices and the complex factor of the exterior polygon mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> prev, C = parameters(m);\n\njulia> prev\n4-element Array{Complex{Float64},1}:\n       1.0+0.0im\n  0.376406-0.926455im\n -0.902383-0.430935im\n -0.186756+0.982406im\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Properties.DerivativeMap-Tuple{}",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Properties.DerivativeMap",
    "category": "method",
    "text": "DerivativeMap(m::ConformalMap)\n\nConstructs new conformal maps from the first and second derivatives of the conformal map m.\n\nThese new conformal maps can be evaluated at a single or vector of points just as  m is. The first entry in the tuple returned is the first derivative, the second entry is the second derivative.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> dm = DerivativeMap(m);\n\njulia> ζ = [0.1,0.5-0.75im,-0.25-0.3im];\n\njulia> dz, ddz = dm(ζ;inside=true);\n\njulia> dz\n3-element Array{Complex{Float64},1}:\n  67.2068+76.6284im\n -1.11666+0.544576im\n  3.99129-5.30641im\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Properties.InverseMap",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Properties.InverseMap",
    "category": "type",
    "text": "InverseMap(m::ConformalMap)\n\nConstructs the inverse conformal map of the conformal map m.\n\nThis inverse conformal map can be evaluated at a single or vector of points. Points should be outside the body. Whether the resulting point in the circle plane is interpreted inside or outside the circle is determined by the optional argument inside, which defaults to false.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> m⁻¹ = InverseMap(m);\n\njulia> ζ = [1.0+3.0im,-2.0-2.0im,0.1+1.1im];\n\njulia> m⁻¹(m(ζ))\n3-element Array{Complex{Float64},1}:\n  1.0+3.0im\n -2.0-2.0im\n  0.1+1.1im\n\n\n\n"
},

{
    "location": "exterior.html#Base.length-Tuple{SchwarzChristoffel.Properties.ConformalMap}",
    "page": "Exterior map",
    "title": "Base.length",
    "category": "method",
    "text": "length(m::ConformalMap) -> Integer\n\nReturns the number of control points/vertices of the map m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> length(m)\n4\n\n\n\n"
},

{
    "location": "exterior.html#Methods-1",
    "page": "Exterior map",
    "title": "Methods",
    "category": "section",
    "text": "Modules = [Exterior]\nOrder   = [:type, :function]"
},

{
    "location": "exterior.html#Index-1",
    "page": "Exterior map",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"exterior.md\"]"
},

]}
