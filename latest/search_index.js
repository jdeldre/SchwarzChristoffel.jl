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
    "text": ""
},

{
    "location": "index.html#About-the-package-1",
    "page": "Home",
    "title": "About the package",
    "category": "section",
    "text": "The purpose of this package is to enable easy construction and evaluation of the mapping from the region inside or outside the unit circle to the exterior of a closed polygon."
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
    "text": "DocTestSetup = quote\nsrand(1)\nendusing SchwarzChristoffelFirst, we create a polygon shape by specifying its vertices. Note that the vertices must be provided in counter-clockwise order.x = [-1.0,0.2,1.0,-1.0]; y = [-1.0,-1.0,0.5,1.0];\np = Polygon(x,y)Let\'s plot the polygon to make sure it matches what we wanted.plot(p)\nsavefig(\"polygon4.svg\",format=\"svg\"); nothing # hide<object data=\"polygon4.svg\" type=\"image/svg+xml\"></object>Now, we create the map from the unit circle to the polygon.m = ExteriorMap(p)Let\'s visualize what we\'ve constructed. Here, we will inspect the mapping from the exterior of the unit circle to the exterior of the polygon.conformal_grid(m)\nsavefig(\"polygongrid.svg\",format=\"svg\"); nothing # hide<object data=\"polygongrid.svg\" type=\"image/svg+xml\"></object>We can now easily evaluate the map at any place we like. It could be evaluated outside the unit circle:zeta = 1.2 + 0.1im\nevaluate(zeta,m)or it could be evaluated inside the unit circle:zeta = 0.5 + 0.1im\nevaluate(zeta,m,true)We can also evaluate the first and second derivative of the map at any place(s). Let\'s evaluate at a range of points outside the circle.zeta = collect(1.1:0.1:2.0) + 0.1im\ndz,ddz = evalderiv(zeta,m)\ndz"
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
    "text": "vertex(p::Polygon) -> -> Vector{Complex128}\n\nReturns the vector of vertices of the polygon p, in complex form.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> vertex(p)\n4-element Array{Complex{Float64},1}:\n -1.0-1.0im\n  0.2-1.0im\n  1.0+0.5im\n -1.0+1.0im\n\n\n\n"
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
    "location": "polygons.html#Methods-1",
    "page": "Polygons",
    "title": "Methods",
    "category": "section",
    "text": "Polygons.Polygon\nPolygons.vertex\nPolygons.interiorangle\nBase.length(::Polygons.Polygon)\nBase.isinf(::Polygons.Polygon)\nPolygons.isinpoly"
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
    "location": "exterior.html#SchwarzChristoffel.Exterior.ExteriorMap",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.ExteriorMap",
    "category": "type",
    "text": "ExteriorMap(p::Polygon[;tol::Float64][,ncoeff::Int])\n\nCreate a Schwarz-Christoffel  from the interior or exterior of the unit circle to the exterior of polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p)\nExterior map with\n   vertices: (-1.0,-1.0), (0.2,-1.0), (1.0,0.5), (-1.0,1.0),\n   interior angles/π: 0.5, 0.656, 0.422, 0.422,\n   prevertices on circle: (1.0,0.0), (0.3764,-0.9265), (-0.9024,-0.4309), (-0.1868,0.9824),\n   prevertex angles/π: -0.7291, -0.3519, 0.1291, 0.7111,\n   constant = 0.6722 + 0.7669im, accuracy = 1.0e-8,\n   number of multipole coefficients = 12\n\nExteriorMap(p;tol=1e-12) manually sets the tolerance to 1e-12 (the default is 1e-8).\n\nExteriorMap(p;ncoeff=16) manually sets the number of coefficients of negative powers of the multipole expansion of the mapping to 16 (the default is 12).\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.parameters",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.parameters",
    "category": "function",
    "text": "parameters(m::ExteriorMap) -> Tuple{Vector{Complex128},Complex128}\n\nReturns a tuple of a vector of the prevertices and the complex factor of the mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> prev, C = parameters(m);\n\njulia> prev\n4-element Array{Complex{Float64},1}:\n       1.0+0.0im\n  0.376406-0.926455im\n -0.902383-0.430935im\n -0.186756+0.982406im\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.coefficients",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.coefficients",
    "category": "function",
    "text": "coefficients(m::ExteriorMap) -> Tuple{Vector{Complex128},Vector{Complex128}}\n\nReturns a tuple of vectors of the complex coefficients of the multipole expansion of the mapping z(zeta) described by m as well as the coefficients of the square magnitude of the mapping z(zeta)^2.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> ccoeff, dcoeff = coefficients(m);\n\njulia> ccoeff\n14-element Array{Complex{Float64},1}:\n       1.0198+0.0im\n    -0.210443-0.015142im\n  -0.00655708+0.0398156im\n     0.136922+0.0951343im\n    -0.095035+0.0891769im\n    0.0184341+0.0299586im\n    0.0136513+2.78095e-5im\n   -0.0159533-0.00264418im\n  -0.00167426-0.00501161im\n  -0.00578705-0.000221652im\n  -0.00447511+0.00252069im\n   0.00469089-0.00150588im\n  0.000441767-0.00192516im\n -0.000381357-0.00174291im \n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.moments",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.moments",
    "category": "function",
    "text": "moments(m::ExteriorMap) -> Vector{Complex128}\n\nReturn the moments of the prevertices for mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> mom = moments(m)\n13-element Array{Complex{Float64},1}:\n -2.46691e-9+3.04899e-9im\n  -0.0128596+0.0780855im\n    0.805587+0.559726im\n    -1.12125+1.04835im\n    0.316471+0.633964im\n    0.462871+0.225702im\n    -1.56266+0.0306815im\n   -0.106975-0.476173im\n   -0.720332-0.0496159im\n     -1.1805+0.0838739im\n      1.1618-0.762023im\n  -0.0612155-0.5728im\n   -0.223423-0.726949im\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.area",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.area",
    "category": "function",
    "text": "area(m::ExteriorMap) -> Float64\n\nReturns the area of the polygon described by the mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> area(m)\n2.9\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.centroid",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.centroid",
    "category": "function",
    "text": "centroid(m::ExteriorMap) -> Complex128\n\nReturns the complex centroid position of the polygon described by the mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> centroid(m)\n-0.20919540229885059 - 0.04022988505747128im\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.Jmoment",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.Jmoment",
    "category": "function",
    "text": "Jmoment(m::ExteriorMap) -> Float64\n\nReturns the second area moment of the polygon described by the mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> Jmoment(m)\n1.5768333333333333\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.evaluate",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.evaluate",
    "category": "function",
    "text": "evaluate(zeta::Vector{Complex128},m::ExteriorMap,inside::Bool) -> Vector{Complex128}\n\nEvaluates the mapping m at the vector of points zeta, which are assumed to lie inside the unit circle if inside is true, or are assumed outside the unit circle if inside is false.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> zeta = [0.1,0.5-0.75im,-0.25-0.3im];\n\njulia> evaluate(zeta,m,true)\n3-element Array{Complex{Float64},1}:\n   -6.9344-7.68965im\n 0.0439774-1.11249im\n   2.41181-0.044779im\n\n\n\nevaluate(zeta::Vector{Complex128},m::ExteriorMap) -> Vector{Complex128}\n\nEvaluates the mapping m at the vector of points zeta, which are assumed to lie outside the unit circle.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];\n\njulia> evaluate(zeta,m)\n3-element Array{Complex{Float64},1}:\n   0.81614+3.02956im\n  -2.25237-2.08523im\n -0.333104+0.975837im\n\n\n\nevaluate(zeta::Complex128,m::ExteriorMap...) -> Complex128\n\nEvaluates m at a single point zeta.\n\n\n\n"
},

{
    "location": "exterior.html#SchwarzChristoffel.Exterior.evalderiv",
    "page": "Exterior map",
    "title": "SchwarzChristoffel.Exterior.evalderiv",
    "category": "function",
    "text": "evalderiv(zeta::Vector{Complex128},m::ExteriorMap,inside::Bool) -> Tuple{Vector{Complex128},Vector{Complex128}}\n\nEvaluates the first and second derivatives of the mapping m at the vector of points zeta, which are assumed to lie inside the unit circle if inside is true, or are assumed outside the unit circle if inside is false. The first entry in the tuple returned is the first derivative, the second entry is the second derivative.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> zeta = [0.1,0.5-0.75im,-0.25-0.3im];\n\njulia> dz, ddz = evalderiv(zeta,m,true);\n\njulia> dz\n3-element Array{Complex{Float64},1}:\n  67.2068+76.6284im\n -1.11666+0.544576im\n  3.99129-5.30641im\n\n\n\nevalderiv(zeta::Vector{Complex128},m::ExteriorMap) -> Tuple{Vector{Complex128},Vector{Complex128}}\n\nEvaluates the first and second derivatives of the mapping m at the vector of points zeta, which are assumed to lie outside the unit circle. The first entry in the tuple returned is the first derivative, the second entry is the second derivative.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];\n\njulia> dz, ddz = evalderiv(zeta,m);\n\njulia> dz\n3-element Array{Complex{Float64},1}:\n 1.03053+0.00444992im\n 1.00696-0.0115011im\n 1.30078-0.266625im\n\n\n\nevalderiv(zeta::Complex128,m::ExteriorMap...) -> Tuple{Complex128,Complex128}\n\nEvaluates the derivatives of m at a single point zeta.\n\n\n\n"
},

{
    "location": "exterior.html#Methods-1",
    "page": "Exterior map",
    "title": "Methods",
    "category": "section",
    "text": "Exterior.ExteriorMap\nExterior.parameters\nExterior.coefficients\nExterior.moments\nExterior.area\nExterior.centroid\nExterior.Jmoment\nExterior.evaluate\nExterior.evalderiv"
},

{
    "location": "exterior.html#Index-1",
    "page": "Exterior map",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"exterior.md\"]"
},

{
    "location": "plotting.html#",
    "page": "Plotting",
    "title": "Plotting",
    "category": "page",
    "text": ""
},

{
    "location": "plotting.html#Plotting-1",
    "page": "Plotting",
    "title": "Plotting",
    "category": "section",
    "text": "DocTestSetup = quote\nusing SchwarzChristoffel\nsrand(1)\nend"
},

{
    "location": "plotting.html#SchwarzChristoffel.plot-Tuple{SchwarzChristoffel.Polygons.Polygon}",
    "page": "Plotting",
    "title": "SchwarzChristoffel.plot",
    "category": "method",
    "text": "plot(p::Polygon)\n\nPlots the polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> plot(p);\n\n\n\n"
},

{
    "location": "plotting.html#SchwarzChristoffel.conformal_grid",
    "page": "Plotting",
    "title": "SchwarzChristoffel.conformal_grid",
    "category": "function",
    "text": "conformal_grid(m::ExteriorMap)\n\nPlots the grid lines generated by the exterior mapping m.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> m = ExteriorMap(p);\n\njulia> conformal_grid(m);\n\n\n\n"
},

{
    "location": "plotting.html#Methods-1",
    "page": "Plotting",
    "title": "Methods",
    "category": "section",
    "text": "plot(::Polygon)\nconformal_grid"
},

{
    "location": "plotting.html#Index-1",
    "page": "Plotting",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"plotting.md\"]"
},

]}
