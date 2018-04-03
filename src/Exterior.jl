module Exterior

using NLsolve
using ..Properties
using ..Polygons
using ..Integration

export ExteriorMap,evaluate,evalderiv,parameters,coefficients,
        moments,area,centroid,Jmoment


struct ExteriorMap <: Map

  vertex :: Vector{Complex128}
  angle  :: Vector{Float64}
  qdata  :: Tuple{Array{Float64,2},Array{Float64,2}}
  accuracy  :: Float64

  prevertex :: Vector{Complex128}
  constant  :: Complex128

  preprev   :: Vector{Complex128}
  prevangle :: Vector{Float64}
  ncoeff    :: Int64
  ccoeff    :: Vector{Complex128}
  dcoeff    :: Vector{Complex128}
  mom       :: Vector{Complex128}

  area      :: Float64
  Zc        :: Complex128
  J         :: Float64

end

"""
    ExteriorMap(p::Polygon[;tol::Float64][,ncoeff::Int])

Create a Schwarz-Christoffel  from the interior or exterior of
the unit circle to the exterior of polygon `p`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p)
Exterior map with
   vertices: (-1.0,-1.0), (0.2,-1.0), (1.0,0.5), (-1.0,1.0),
   interior angles/π: 0.5, 0.656, 0.422, 0.422,
   prevertices on circle: (1.0,0.0), (0.3764,-0.9265), (-0.9024,-0.4309), (-0.1868,0.9824),
   prevertex angles/π: -0.7291, -0.3519, 0.1291, 0.7111,
   constant = 0.6722 + 0.7669im, accuracy = 1.0e-8,
   number of multipole coefficients = 12
```

`ExteriorMap(p;tol=1e-12)` manually sets the tolerance to `1e-12`
(the default is 1e-8).

`ExteriorMap(p;ncoeff=16)` manually sets the number of coefficients of
negative powers of the multipole expansion of the mapping to `16`
(the default is 12).
"""
function ExteriorMap(p::Polygon;tol::Float64 = 1e-8,ncoeff::Int = 12)

  n = length(p.vert)

  zeta0 = Complex128[]

  w = flipdim(vertex(p),1)
  beta = 1.-flipdim(interiorangle(p),1)

  # do some fixing
  n = length(p)
  renum = 1:n
  shift = [2:n;1]
  w = flipdim(vertex(p),1)
  beta = 1 - flipdim(interiorangle(p),1)
  while any(abs.(beta[n]-[0;1]).<eps()) & (n > 2)
    renum = renum[shift]
    w = w[shift]
    beta = beta[shift]
    if renum[1]==1
      deg = abs.(beta) .< eps()
      w[deg] = []
      beta[deg] = []
      renum = 1:2
      n = 2
      break
    end
  end

  p = Polygon(flipdim(w,1),1.-flipdim(beta,1))

  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)

  zeta, c = param(w,beta,zeta0,qdat)
  preprev = -c/abs(c)./flipdim(zeta,1)
  prevangle = angle.(preprev)

  # compute geometric stuff
  z = vertex(p)
  zmid = 0.5*(z+circshift(z,-1))
  dz = circshift(z,-1)-z
  area = 0.5*imag(sum(conj.(zmid).*dz))
  if area > 0
    Zc = -0.5im/area*sum(dz.*(abs.(zmid).^2+abs.(dz).^2/12))
  else
    Zc = mean(z)
  end
  J = 0.25*imag(sum(conj.(zmid).*dz.*(abs.(zmid).^2+abs.(dz).^2/12)))

  # first two entries are for c₁ and c₀.
  betaflip = flipdim(beta,1)
  ccoeff = Complex128[abs(c),0.0]
  mom = [sum(betaflip.*preprev)]
  for k = 1:ncoeff
    push!(mom,sum(betaflip.*preprev.^(k+1)))
    coeffk = abs(c)*getcoefflist(k+1,1,mom);
    push!(ccoeff,coeffk)
  end

  # fix a
  zetainf = 5.0*exp.(im*collect(0:π/2:3π/2))
  sigmainf = -c/abs(c)./zetainf
  zinf = evaluate(sigmainf,w,beta,zeta,c,qdat)
  ccoeff[2] = mean([zinf[i]-sum(ccoeff.*zetainf[i].^(1:-1:-ncoeff))
                  for i = 1:length(zetainf)])

  # first entry is d₀
  dcoeff = [dot(ccoeff,ccoeff)]
  for k = 1:ncoeff+1
    push!(dcoeff,dot(ccoeff[1:end-k],ccoeff[k+1:end]))
  end

  ExteriorMap(p.vert,p.angle,qdat,tol,zeta,c,preprev,prevangle,
              ncoeff,ccoeff,dcoeff,mom,area,Zc,J)
end

function Base.show(io::IO, m::ExteriorMap)
    println(io, "Exterior map with")
    print(io,   "   ")
    print(io,"vertices: ")
    for i = 1:length(m.vertex)
        print(io,"($(round(real(m.vertex[i]),4)),$(round(imag(m.vertex[i]),4))), ")
    end
    println(io)
    #for i = 1:length(m.vertex)
    #println(io, "       $(round(m.vertex[i],4))")
    #end
    print(io,"   ")
    print(io, "interior angles/π: ")
    for i = 1:length(m.angle)
        print(io, "$(round(m.angle[i],4)), ")
    end
    println(io)
    print(io,"   ")
    print(io, "prevertices on circle: ")
    for i = length(m.prevertex):-1:1
        print(io,"($(round(real(m.prevertex[i]),4)),$(round(imag(m.prevertex[i]),4))), ")
    end
    println(io)
    print(io, "   ")
    print(io, "prevertex angles/π: ")
    for i = 1:length(m.prevangle)
        print(io, "$(round(m.prevangle[i]/π,4)), ")
    end
    println(io)
    print(io, "   ")
    print(io, "constant = $(round(m.constant,4)), ")
    print(io, "accuracy = $(m.accuracy), ")
    println(io)
    print(io, "   ")
    print(io, "number of multipole coefficients = $(m.ncoeff)")
    println(io)

end

"""
    parameters(m::ExteriorMap) -> Tuple{Vector{Complex128},Complex128}

Returns a tuple of a vector of the prevertices and the complex factor of
the mapping `m`.

# Example

```jldoctest paramtest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> prev, C = parameters(m);

julia> prev
4-element Array{Complex{Float64},1}:
       1.0+0.0im
  0.376406-0.926455im
 -0.902383-0.430935im
 -0.186756+0.982406im
```
"""


parameters(m::ExteriorMap) = flipdim(m.prevertex,1), m.constant

"""
    coefficients(m::ExteriorMap) -> Tuple{Vector{Complex128},Vector{Complex128}}

Returns a tuple of vectors of the complex coefficients of the multipole
expansion of the mapping \$z(\\zeta)\$ described by `m` as well as the
coefficients of the square magnitude of the mapping \$|z(\\zeta)|^2\$.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> ccoeff, dcoeff = coefficients(m);

julia> ccoeff
14-element Array{Complex{Float64},1}:
       1.0198+0.0im
    -0.210364-0.0161983im
  -0.00655708+0.0398156im
     0.136922+0.0951343im
    -0.095035+0.0891769im
    0.0184341+0.0299586im
    0.0136513+2.78095e-5im
   -0.0159533-0.00264418im
  -0.00167426-0.00501161im
  -0.00578705-0.000221652im
  -0.00447511+0.00252069im
   0.00469089-0.00150588im
  0.000441767-0.00192516im
 -0.000381357-0.00174291im
```
"""
coefficients(m::ExteriorMap) = m.ccoeff, m.dcoeff

"""
    moments(m::ExteriorMap) -> Vector{Complex128}

Return the moments of the prevertices for mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> mom = moments(m)
13-element Array{Complex{Float64},1}:
 -2.46691e-9+3.04899e-9im
  -0.0128596+0.0780855im
    0.805587+0.559726im
    -1.12125+1.04835im
    0.316471+0.633964im
    0.462871+0.225702im
    -1.56266+0.0306815im
   -0.106975-0.476173im
   -0.720332-0.0496159im
     -1.1805+0.0838739im
      1.1618-0.762023im
  -0.0612155-0.5728im
   -0.223423-0.726949im
```
"""
moments(m::ExteriorMap) = m.mom

"""
    area(m::ExteriorMap) -> Float64

Returns the area of the polygon described by the mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> area(m)
2.9
```
"""
area(m::ExteriorMap) = m.area

"""
    centroid(m::ExteriorMap) -> Complex128

Returns the complex centroid position of the polygon described by the
mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> centroid(m)
-0.20919540229885059 - 0.04022988505747128im
```
"""
centroid(m::ExteriorMap) = m.Zc

"""
    Jmoment(m::ExteriorMap) -> Float64

Returns the second area moment of the polygon described by the
mapping `m`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> Jmoment(m)
1.5768333333333333
```
"""
Jmoment(m::ExteriorMap) = m.J


function param(w::Vector{Complex128},beta::Vector{Float64},
                 zeta0::Vector{Complex128},
                 qdat::Tuple{Array{Float64,2},Array{Float64,2}})
  # w clockwise
  # beta turning angles

  n = length(w)
  if n == 2
    zeta = Complex128[-1,1]
  else
    len = abs.(diff(circshift(w,1)))
    nmlen = abs.(len[3:n-1]/len[2])
    if isempty(zeta0)
      y0 = zeros(n-1)
    else
      zeta0 = zeta0/zeta0[n]
      θ = angle.(zeta0)
      θ[θ.<=0] = θ[θ.<=0] + 2π
      dt = diff([0;θ[1:n-1];2π])
      @. y0 = log(dt[1:n-1]/dt[2:n])
    end

    depfun!(F,y) = depfunfull!(F,y,n,beta,nmlen,qdat)

    F0 = similar(y0)
    df = OnceDifferentiable(depfun!, y0, F0)
    sol = nlsolve(depfun!,y0,show_trace = :false)

    zeta, θ = y_to_zeta(sol.zero)

  end

  mid = zeta[1]*exp(0.5*im*angle(zeta[2]/zeta[1]))
  dequad = DQuad(zeta,beta,qdat)
  c = (w[2]-w[1])/(dequad(zeta[1],mid,1)-dequad(zeta[2],mid,2))

  return zeta, c

end

function evaluate(zeta::Vector{Complex128},w::Vector{Complex128},
                  beta::Vector{Float64},prev::Vector{Complex128},
                  c::Complex128,qdat::Tuple{Array{Float64,2},Array{Float64,2}})

  # this assumes zeta inside the unit circle

  if isempty(zeta)
    nothing
  end

  n = length(w)
  neval = length(zeta)
  tol = 10.0^(-size(qdat[1],1))

  # set up the integrator
  dequad = DQuad(prev,beta,qdat)

  # initialize the mapped evaluation points
  wp = zeros(Complex128,neval)

  # find the closest prevertices to each evaluation point and their
  #  corresponding distances
  dz = abs.(hcat([zeta for i=1:n]...)-vcat([transpose(prev) for i=1:neval]...))
  (dist,ind) = findmin(dz,2)
  sing = floor.(Int,(ind[:]-1)/neval)+1

  # find any prevertices in the evaluation list and set them equal to
  #  the corresponding vertices. The origin is also a singular point
  vertex = (dist[:] .< tol)
  wp[vertex] = w[sing[vertex]]
  zerop = abs.(zeta) .< tol
  wp[zerop] = Inf
  vertex = vertex .| zerop

  # the starting (closest) singularities for each evaluation point
  prevs = prev[sing]

  # set the initial values of the non-vertices
  wp[.!vertex] = w[sing[.!vertex]]

  # distance to singularity at the origin
  abszeta = abs.(zeta)

  # unfinished cases
  unf = .!vertex

  # set the integration starting points
  zetaold = copy(prevs)
  zetanew = copy(zetaold)
  dist = ones(neval)
  while any(unf)
    # distance to integrate still
    dist[unf] = min.(1,2*abszeta[unf]./abs.(zeta[unf]-zetaold[unf]))

    # new integration end point
    zetanew[unf] = zetaold[unf] + dist[unf].*(zeta[unf]-zetaold[unf])

    # integrate
    wp[unf] = wp[unf] + c*dequad.(zetaold[unf],zetanew[unf],sing[unf])

    # set new starting integration points for those that can be integrated
    #  further
    unf = dist .< 1
    zetaold[unf] = zetanew[unf]

    # only the first step can have a singularity
    sing .= 0

  end

  return wp

end

function evalderiv(zeta::Vector{Complex128},beta::Vector{Float64},
                  prev::Vector{Complex128},c::Complex128)

    n = length(prev)
    neval = length(zeta)
    beta = [beta;-2]
    terms = [hcat([1 .- zeta/prev[i] for i = 1:n]...) zeta]
    dz = c*exp.(log.(terms)*beta)
    terms2 = [[-beta[i]/prev[i] for i = 1:n];beta[n+1]]
    ddz = dz.*((1.0./terms)*terms2)
    return dz, ddz

end

"""
    evaluate(zeta::Vector{Complex128},m::ExteriorMap,inside::Bool) -> Vector{Complex128}

Evaluates the mapping `m` at the vector of points `zeta`, which are
assumed to lie inside the unit circle if `inside` is `true`, or
are assumed outside the unit circle if `inside` is `false`.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> zeta = [0.1,0.5-0.75im,-0.25-0.3im];

julia> evaluate(zeta,m,true)
3-element Array{Complex{Float64},1}:
   -6.9344-7.68965im
 0.0439774-1.11249im
   2.41181-0.044779im
```
"""
function evaluate(zeta::Vector{Complex128},m::ExteriorMap,inside::Bool)

  if inside
    return evaluate(zeta,flipdim(m.vertex,1),1.-flipdim(m.angle,1),
            m.prevertex,m.constant,m.qdata)
  else
    b = -m.constant/abs(m.constant)
    zeta[zeta.==0] = eps();
    zeta[abs.(zeta).<1] = zeta[abs.(zeta).<1]./abs.(zeta[abs.(zeta).<1])

    sigma = b./zeta
    return evaluate(sigma,flipdim(m.vertex,1),1.-flipdim(m.angle,1),
            m.prevertex,m.constant,m.qdata)
  end
end

"""
    evaluate(zeta::Vector{Complex128},m::ExteriorMap) -> Vector{Complex128}

Evaluates the mapping `m` at the vector of points `zeta`, which are
assumed to lie outside the unit circle.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];

julia> evaluate(zeta,m)
3-element Array{Complex{Float64},1}:
   0.81614+3.02956im
  -2.25237-2.08523im
 -0.333104+0.975837im
```
"""
evaluate(zeta::Vector{Complex128},m::ExteriorMap) = evaluate(zeta,m,false)

"""
    evaluate(zeta::Complex128,m...) -> Complex128

Evaluates `m` at a single point `zeta`.
"""
evaluate(zeta::Complex128,m...) = getindex(evaluate([zeta],m...),1)

"""
    evalderiv(zeta::Vector{Complex128},m::ExteriorMap,inside::Bool) -> Tuple{Vector{Complex128},Vector{Complex128}}

Evaluates the first and second derivatives of the mapping `m` at the vector
of points `zeta`, which are assumed to lie inside the unit circle if
`inside` is `true`, or are assumed outside the unit circle if `inside` is
`false`. The first entry in the tuple returned is the first derivative,
the second entry is the second derivative.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> zeta = [0.1,0.5-0.75im,-0.25-0.3im];

julia> dz, ddz = evalderiv(zeta,m,true);

julia> dz
3-element Array{Complex{Float64},1}:
  67.2068+76.6284im
 -1.11666+0.544576im
  3.99129-5.30641im
```
"""
function evalderiv(zeta::Vector{Complex128},m::ExteriorMap,inside::Bool)

  if inside
    return evalderiv(zeta,1.-flipdim(m.angle,1),m.prevertex,m.constant)
  else
    b = -m.constant/abs(m.constant)
    zeta[zeta.==0] = eps();
    zeta[abs.(zeta).<1] = zeta[abs.(zeta).<1]./abs.(zeta[abs.(zeta).<1])

    sigma = b./zeta
    dsigma = -sigma./zeta
    ddsigma = -2.0*dsigma./zeta
    dz, ddz = evalderiv(sigma,1.-flipdim(m.angle,1),m.prevertex,m.constant)
    ddz = ddz.*dsigma.^2 + dz.*ddsigma
    dz .*= dsigma
    return dz, ddz
  end
end

"""
    evalderiv(zeta::Vector{Complex128},m::ExteriorMap) -> Tuple{Vector{Complex128},Vector{Complex128}}

Evaluates the first and second derivatives of the mapping `m` at the vector
of points `zeta`, which are assumed to lie outside the unit circle.
The first entry in the tuple returned is the first derivative,
the second entry is the second derivative.

# Example

```jldoctest
julia> p = Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0]);

julia> m = ExteriorMap(p);

julia> zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im];

julia> dz, ddz = evalderiv(zeta,m);

julia> dz
3-element Array{Complex{Float64},1}:
 1.03053+0.00444992im
 1.00696-0.0115011im
 1.30078-0.266625im
```
"""
evalderiv(zeta::Vector{Complex128},m::ExteriorMap) = evalderiv(zeta,m,false)

"""
    evalderiv(zeta::Complex128,m...) -> Tuple{Complex128,Complex128}

Evaluates the derivatives of `m` at a single point `zeta`.
"""
evalderiv(zeta::Complex128,m...) =
      getindex(getindex(evalderiv([zeta],m...),1),1),
          getindex(getindex(evalderiv([zeta],m...),2),1)



function getcoefflist(power,div,mom)

  # mom are the moments, mom[1] is M₁, etc

  if power%div!=0
    error("Indivisible power")
  end
  pow = Int(power/div)

  # Find the set of multi-indices I for which
  # sum(k(t)*t) = power/div. Each row of I corresponds
  # to a different multi-index in the set
  I = [1]
  for j = 2:pow
    I = madvance(I)
  end

  # Find the coefficient for 1/zeta^(pow-1) in the fhat expansion
  coeff = 0
  for j = 1:size(I,1)
    sumI = sum(I[j,:])
    fact = 1
    for l = 1:pow
      il = I[j,l]
      fact *= mom[l]^il/factorial(il)/l^il
    end
    coeff += fact*(-1)^sumI
  end
  return -coeff/(pow-1)

end


function madvance(P)
  # Given a set P of multi-indices
  # that satisfy a certain moment condition
  # |p|_1 = m (for all p in P), find the set Pplus1
  # that satisfies |p'|_1 = m+1 (for all p' in Pplus1)
  nP = size(P,1)
  m = size(P,2)
  Ppad = [ones(nP) P zeros(nP)]
  rows = 0
  Pplus1 = Int64[]
  for k = 1:nP
    # Loop through each index in the kth
    #  member of P, and find non-zero indices for shifting
    for j = m:-1:0
      if Ppad[k,j+1] == 1
        rows += 1
        Pplus1 = [Pplus1;zeros(Int64,1,m+2)]
        Pplus1[rows,1:m+2] = Ppad[k,1:m+2]
        Pplus1[rows,j+2] = Pplus1[rows,j+2]+1
        Pplus1[rows,j+1] = Pplus1[rows,j+1]-1
      end
    end
  end
  Pplus1 = Pplus1[:,2:end]
  Pplus1 = sortrows(Pplus1)
  dP = sum(abs.([transpose(Pplus1[1,:]);diff(Pplus1,1)]),2)
  return Pplus1[dP[:].!=0,:]

end

struct DabsQuad{n}
  zeta :: Vector{Complex128}
  beta :: Vector{Float64}
  nqpts :: Int64
  qdat :: Tuple{Array{Float64,2},Array{Float64,2}}
end

function DabsQuad(zeta::Vector{Complex128},beta::Vector{Float64},tol::Float64)
  n = length(zeta)
  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)
  DabsQuad{n}(zeta,beta,nqpts,qdat)
end

function DabsQuad(zeta::Vector{Complex128},beta::Vector{Float64},
                  qdat::Tuple{Array{Float64,2},Array{Float64,2}})
  n = length(zeta)
  nqpts = size(qdat[1],1)
  DabsQuad{n}(zeta,beta,nqpts,qdat)
end

struct DQuad{n}
  zeta :: Vector{Complex128}
  beta :: Vector{Float64}
  nqpts :: Int64
  qdat :: Tuple{Array{Float64,2},Array{Float64,2}}
end

function DQuad(zeta::Vector{Complex128},beta::Vector{Float64},tol::Float64)
  n = length(zeta)
  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)
  DQuad{n}(zeta,beta,nqpts,qdat)
end

function DQuad(zeta::Vector{Complex128},beta::Vector{Float64},
               qdat::Tuple{Array{Float64,2},Array{Float64,2}})
  n = length(zeta)
  nqpts = size(qdat[1],1)
  DQuad{n}(zeta,beta,nqpts,qdat)
end


function (I::DabsQuad{n})(zeta1::Complex128,zeta2::Complex128,sing1::Int64) where n

   (qnode,qwght) = I.qdat
   argz = angle.(I.zeta)

   argz1 = angle(zeta1)
   argz2 = angle(zeta2)
   ang21 = angle(zeta2/zeta1)

   if (argz2-argz1)*ang21 < 0
      argz2 += 2π*sign(ang21)
   end
   if isempty(sing1)
     sing1 = 0
   end
   result = 0.0
   if zeta1 != zeta2
     zetas = [I.zeta[1:sing1-1];I.zeta[sing1+1:end]]
     dist = min(1,2*minimum(abs.(zetas-zeta1))/abs(zeta2-zeta1))
     argr = argz1 + dist*(argz2-argz1)
     ind = ((sing1+n) % (n+1)) + 1
     nd = 0.5*((argr-argz1)*qnode[:,ind] + argr + argz1)
     wt = 0.5*abs(argr-argz1)*qwght[:,ind]
     θ = hcat([(nd - argz[i] + 2π).%(2π) for i = 1:n]...)
     θ[θ.>π] = 2π-θ[θ.>π]
     terms = 2sin.(0.5θ)
     if !any(terms==0.0)
        if sing1 > 0
            terms[:,sing1] ./= abs.(nd-argz1)
            wt .*= (0.5*abs.(argr-argz1)).^I.beta[sing1]
        end
        result = transpose(exp.(log.(terms)*I.beta))*wt
        while dist < 1.0
            argl = argr
            zetal = exp(im*argl)
            dist = min(1,2*minimum(abs.(I.zeta-zetal)/abs(zetal-zeta2)))
            argr = argl + dist*(argz2-argl)
            nd = 0.5*((argr-argl)*qnode[:,n+1] + argr + argl)
            wt = 0.5*abs(argr-argl)*qwght[:,n+1]
            θ = hcat([(nd - argz[i] + 2π).%(2π) for i = 1:n]...)
            θ[θ.>π] = 2π-θ[θ.>π]
            terms = 2sin.(0.5θ)
            result += transpose(exp.(log.(terms)*I.beta))*wt
        end
    end
   end
   return result
end

function (I::DQuad{n})(zeta1::Complex128,zeta2::Complex128,sing1::Int64) where n

   (qnode,qwght) = I.qdat

   beta = [I.beta;-2]

   if isempty(sing1)
     sing1 = 0
   end
   result = Complex128[0]
   if zeta1 != zeta2
     zetas = [I.zeta[1:sing1-1];I.zeta[sing1+1:end]]
     dist = min(1,2*minimum(abs.(zetas-zeta1))/abs(zeta2-zeta1))
     zetar = zeta1 + dist*(zeta2-zeta1)
     ind = sing1 + (n+1)*(sing1==0)

     nd = 0.5*((zetar-zeta1)*qnode[:,ind] + zetar + zeta1)
     wt = 0.5*(zetar-zeta1)*qwght[:,ind]
     terms = hcat([1 .- nd/I.zeta[i] for i = 1:n]...)
     if !any(terms==0.0)
       terms = hcat(terms,nd)
        if sing1 > 0
            terms[:,sing1] ./= abs.(terms[:,sing1])
            wt .*= (0.5*abs.(zetar-zeta1)).^beta[sing1]
        end
        result = transpose(exp.(log.(terms)*beta))*wt
        while dist < 1.0
            zetal = zetar
            dist = min(1,2*minimum(abs.(I.zeta-zetal)/abs(zetal-zeta2)))
            zetar = zetal + dist*(zeta2-zetal)
            nd = 0.5*((zetar-zetal)*qnode[:,n+1] + zetar + zetal)
            wt = 0.5*(zetar-zetal)*qwght[:,n+1]
            terms = hcat([1 .- nd/I.zeta[i] for i = 1:n]...)
            terms = hcat(terms,nd)
            result += transpose(exp.(log.(terms)*beta))*wt
        end
      end
    end
   return result
end




function depfunfull!(F,y,n,beta,nmlen,qdat)

  zeta, θ = y_to_zeta(y)

  mid = exp.(im*0.5*(θ[1:n-2]+θ[2:n-1]))

  dabsquad = DabsQuad(zeta,beta,qdat)
  #dabsquad(z1,z2,sing1) = dabsquad(z1,z2,sing1,zeta,beta,qdat)

  ints = dabsquad.(zeta[1:n-2],mid,collect(1:n-2))+dabsquad.(zeta[2:n-1],mid,collect(2:n-1))

  if n > 3
    F[1:n-3] = abs.(ints[2:n-2])/abs(ints[1]) - nmlen
  end

  res = -sum(beta./zeta)/ints[1]
  F[n-2] = real(res)
  F[n-1] = imag(res)

end

function y_to_zeta(y::Vector{Float64})

  cs = cumsum(cumprod([1;exp.(-y)]))
  n = length(cs)
  θ = 2π*cs[1:n-1]./cs[n]
  zeta = ones(Complex128,n)
  zeta[1:n-1] = exp.(im*θ)
  return zeta, θ

end

end
