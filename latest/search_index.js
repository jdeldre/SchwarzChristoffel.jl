var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "SchwarzChristoffel",
    "title": "SchwarzChristoffel",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#SchwarzChristoffel-1",
    "page": "SchwarzChristoffel",
    "title": "SchwarzChristoffel",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#About-the-package-1",
    "page": "SchwarzChristoffel",
    "title": "About the package",
    "category": "section",
    "text": "The purpose of this package is to enable easy construction and evaluation of the mapping from the region inside or outside the unit circle to the exterior of a closed polygon."
},

{
    "location": "index.html#Installation-1",
    "page": "SchwarzChristoffel",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia 0.6- and above. It is not a registered package, so it should be installed with:julia> Pkg.clone(\"git@github.com:jdeldre/SchwarzChristoffel.jl.git\")Since it is still under heavy development, you should runjulia> Pkg.test(\"SchwarzChristoffel\") # might take some timeto make sure things are working as intended andjulia> Pkg.update()to get the most recent version of the library and its dependencies.Examples can be found in the documentation and the Jupyter notebooks."
},

{
    "location": "reference.html#",
    "page": "Reference page",
    "title": "Reference page",
    "category": "page",
    "text": ""
},

{
    "location": "reference.html#Reference-page-1",
    "page": "Reference page",
    "title": "Reference page",
    "category": "section",
    "text": "DocTestSetup = quote\nusing SchwarzChristoffel\nsrand(1)\nend"
},

{
    "location": "reference.html#SchwarzChristoffel.Polygons.Polygon",
    "page": "Reference page",
    "title": "SchwarzChristoffel.Polygons.Polygon",
    "category": "type",
    "text": "Polygons.Polygon(x::Vector{Float64}, y::Vector{Float64})\n\nSets up a polygon with the coordinates of the vertices specified with vectors x and y.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\n\n\nPolygons.Polygon(w::Vector{Complex128})\n\nSets up a polygon with the coordinates of the vertices specified with complex vector w.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0-1.0im,0.2-1.0im,1.0+0.5im,-1.0+1.0im])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\n\n\n"
},

{
    "location": "reference.html#SchwarzChristoffel.Polygons.vertex",
    "page": "Reference page",
    "title": "SchwarzChristoffel.Polygons.vertex",
    "category": "function",
    "text": "Polygons.vertex(p::Polygons.Polygon)\n\nReturns the vector of vertices of the polygon p, in complex form.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\njulia> Polygons.vertex(p)\n4-element Array{Complex{Float64},1}:\n -1.0-1.0im\n  0.2-1.0im\n  1.0+0.5im\n -1.0+1.0im\n\n\n\n"
},

{
    "location": "reference.html#SchwarzChristoffel.Polygons.interiorangle",
    "page": "Reference page",
    "title": "SchwarzChristoffel.Polygons.interiorangle",
    "category": "function",
    "text": "Polygons.interiorangle(p::Polygons.Polygon)\n\nReturns the vector of interior angles of the polygon p.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\njulia> Polygons.interiorangle(p)\n4-element Array{Float64,1}:\n 0.5\n 0.655958\n 0.422021\n 0.422021\n\n\n\n"
},

{
    "location": "reference.html#Base.length-Tuple{SchwarzChristoffel.Polygons.Polygon}",
    "page": "Reference page",
    "title": "Base.length",
    "category": "method",
    "text": "length(p::Polygons.Polygon)\n\nReturns the number of vertices of the polygon p.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\njulia> length(p)\n4\n\n\n\n"
},

{
    "location": "reference.html#Base.isinf-Tuple{SchwarzChristoffel.Polygons.Polygon}",
    "page": "Reference page",
    "title": "Base.isinf",
    "category": "method",
    "text": "isinf(p::Polygons.Polygon)\n\nReturns true if any vertex in polygon p is at infinity.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\njulia> isinf(p)\nfalse\n\n\n\n"
},

{
    "location": "reference.html#SchwarzChristoffel.Polygons.isinpoly",
    "page": "Reference page",
    "title": "SchwarzChristoffel.Polygons.isinpoly",
    "category": "function",
    "text": "Polygons.isinpoly(z::Complex128,p::Polygons.Polygon)\n\nReturns true or false depending on whether z is inside or outside polygon p.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\njulia> Polygons.isinpoly(0.0+0.0im,p)\ntrue\n\njulia> Polygons.isinpoly(1.0+2.0im,p)\nfalse\n\nPolygons.isinpoly(z::Complex128,p::Polygons.Polygon,tol::Float64)\n\nReturns true if z is inside or within distance tol polygon p.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\njulia> Polygons.isinpoly(-1.01+0.0im,p)\nfalse\n\njulia> Polygons.isinpoly(-1.01+0.0im,p,1e-2)\ntrue\n\n\n\n"
},

{
    "location": "reference.html#SchwarzChristoffel.Polygons.plot",
    "page": "Reference page",
    "title": "SchwarzChristoffel.Polygons.plot",
    "category": "function",
    "text": "Polygons.plot(p::Polygons.Polygon)\n\nPlots the polygon p.\n\nExample\n\njulia> p = Polygons.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])\nPolygon with 4 vertices at Complex{Float64}[-1.0-1.0im, 0.2-1.0im, 1.0+0.5im, -1.0+1.0im]\n             interior angles/π = [0.5, 0.656, 0.422, 0.422]\n\njulia> Polygons.plot(p);\n\n\n\n"
},

{
    "location": "reference.html#Methods-1",
    "page": "Reference page",
    "title": "Methods",
    "category": "section",
    "text": "Polygons.Polygon\nPolygons.vertex\nPolygons.interiorangle\nBase.length(::Polygons.Polygon)\nBase.isinf(::Polygons.Polygon)\nPolygons.isinpoly\nPolygons.plot"
},

{
    "location": "reference.html#Index-1",
    "page": "Reference page",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"reference.md\"]"
},

]}
