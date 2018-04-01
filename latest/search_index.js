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
    "text": "Polygon(x::Vector{Float64}, y::Vector{Float64})\n\nA polygon defined by its vertices and associated interior angles.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at\n             (-1.0,-1.0) (0.2,-1.0) (1.0,0.5) (-1.0,1.0)\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\n\n\n"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.vertex",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.vertex",
    "category": "function",
    "text": "vertex(p::Polygon)\n\nReturns the vector of vertices of the polygon p, in complex form.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> vertex(p)\n4-element Array{Complex{Float64},1}:\n -1.0-1.0im\n  0.2-1.0im\n  1.0+0.5im\n -1.0+1.0im\n\n\n\n"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.interiorangle",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.interiorangle",
    "category": "function",
    "text": "interiorangle(p::Polygon)\n\nReturns the vector of interior angles (divided by pi) of the polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> interiorangle(p)\n4-element Array{Float64,1}:\n 0.5\n 0.655958\n 0.422021\n 0.422021\n\n\n\n"
},

{
    "location": "polygons.html#Base.length-Tuple{SchwarzChristoffel.Polygons.Polygon}",
    "page": "Polygons",
    "title": "Base.length",
    "category": "method",
    "text": "length(p::Polygon)\n\nReturns the number of vertices of the polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> length(p)\n4\n\n\n\n"
},

{
    "location": "polygons.html#Base.isinf-Tuple{SchwarzChristoffel.Polygons.Polygon}",
    "page": "Polygons",
    "title": "Base.isinf",
    "category": "method",
    "text": "isinf(p::Polygon)\n\nReturns true if any vertex in polygon p is at infinity.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> isinf(p)\nfalse\n\n\n\n"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.isinpoly",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.isinpoly",
    "category": "function",
    "text": "isinpoly(z::Complex128,p::Polygon)\n\nReturns true or false depending on whether z is inside or outside polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> isinpoly(0.0+0.0im,p)\ntrue\n\njulia> isinpoly(1.0+2.0im,p)\nfalse\n\n\n\nisinpoly(z::Complex128,p::Polygon,tol::Float64)\n\nReturns true if z is inside or within distance tol of polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> isinpoly(-1.01+0.0im,p)\nfalse\n\njulia> isinpoly(-1.01+0.0im,p,1e-2)\ntrue\n\n\n\n"
},

{
    "location": "polygons.html#SchwarzChristoffel.Polygons.plot",
    "page": "Polygons",
    "title": "SchwarzChristoffel.Polygons.plot",
    "category": "function",
    "text": "plot(p::Polygon)\n\nPlots the polygon p.\n\nExample\n\njulia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);\n\njulia> plot(p);\n\n\n\n"
},

{
    "location": "polygons.html#Methods-1",
    "page": "Polygons",
    "title": "Methods",
    "category": "section",
    "text": "Polygons.Polygon\nPolygons.vertex\nPolygons.interiorangle\nBase.length(::Polygons.Polygon)\nBase.isinf(::Polygons.Polygon)\nPolygons.isinpoly\nPolygons.plot"
},

{
    "location": "polygons.html#Index-1",
    "page": "Polygons",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"polygons.md\"]"
},

]}
