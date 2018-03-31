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

function ExteriorMap(p::Polygon;tol::Float64 = 1e-8,ncoeff::Int = 12)

  n = length(p.vert)

  zeta0 = Complex128[]

  w = flipdim(vertex(p),1)
  beta = 1.-flipdim(interiorangle(p),1)

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
    Zc = mean.(z)
  end
  J = 0.25*imag(sum(conj.(zmid).*dz.*(abs.(zmid).^2+abs.(dz).^2/12)))

  zetainf = Complex128[100.0]
  sigmainf = -c/abs(c)./zetainf
  zinf = evaluate(sigmainf,w,beta,zeta,c,qdat)
  a = zinf[1] - abs(c)*zetainf[1]

  # first two entries are for c₁ and c₀.
  beta = flipdim(beta,1)
  ccoeff = Complex128[abs(c),a]
  mom = [sum(beta.*preprev)]
  for k = 1:ncoeff
    push!(mom,sum(beta.*preprev.^(k+1)))
    coeffk = abs(c)*getcoefflist(k+1,1,mom);
    push!(ccoeff,coeffk)
  end

  # first entry is d₀
  dcoeff = [dot(ccoeff,ccoeff)]
  for k = 1:ncoeff+1
    push!(dcoeff,dot(ccoeff[1:end-k],ccoeff[k+1:end]))
  end

  ExteriorMap(p.vert,p.angle,qdat,tol,zeta,c,preprev,prevangle,
              ncoeff,ccoeff,dcoeff,mom,area,Zc,J)
end

function Base.show(io::IO, map::ExteriorMap)
    println(io, "Exterior map with $(length(map.vertex)) vertices at")
    for i = 1:length(map.vertex)
    println(io, "       $(round(map.vertex[i],4))")
    end
    println(io, "   interior angles/π")
    for i = 1:length(map.angle)
    println(io, "       $(round(map.angle[i],4))")
    end
    println(io, "   prevertices on circle at")
    for i = length(map.prevertex):-1:1
    println(io, "       $(round(map.prevertex[i],4))")
    end
    println(io, "   constant = $(round(map.constant,4))")
    println(io, "   accuracy = $(map.accuracy)")
    println(io, "   number of multipole coefficients = $(map.ncoeff)")

end

parameters(map::ExteriorMap) = flipdim(map.prevertex,1), map.constant

coefficients(map::ExteriorMap) = map.ccoeff, map.dcoeff

moments(map::ExteriorMap) = map.mom

area(map::ExteriorMap) = map.area
centroid(map::ExteriorMap) = map.Zc
Jmoment(map::ExteriorMap) = map.J


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

function evaluate(zeta::Vector{Complex128},map::ExteriorMap,inside::Bool)

  if inside
    return evaluate(zeta,flipdim(map.vertex,1),1.-flipdim(map.angle,1),
            map.prevertex,map.constant,map.qdata)
  else
    b = -map.constant/abs(map.constant)
    zeta[zeta.==0] = eps();
    zeta[abs.(zeta).<1] = zeta[abs.(zeta).<1]./abs.(zeta[abs.(zeta).<1])

    sigma = b./zeta
    return evaluate(sigma,flipdim(map.vertex,1),1.-flipdim(map.angle,1),
            map.prevertex,map.constant,map.qdata)
  end
end

evaluate(zeta::Vector{Complex128},map::ExteriorMap) = evaluate(zeta,map,false)

function evalderiv(zeta::Vector{Complex128},map::ExteriorMap,inside::Bool)

  if inside
    return evalderiv(zeta,1.-flipdim(map.angle,1),map.prevertex,map.constant)
  else
    b = -map.constant/abs(map.constant)
    zeta[zeta.==0] = eps();
    zeta[abs.(zeta).<1] = zeta[abs.(zeta).<1]./abs.(zeta[abs.(zeta).<1])

    sigma = b./zeta
    dsigma = -sigma./zeta
    ddsigma = -2.0*dsigma./zeta
    dz, ddz = evalderiv(sigma,1.-flipdim(map.angle,1),map.prevertex,map.constant)
    ddz = ddz.*dsigma.^2 + dz.*ddsigma
    dz .*= dsigma
    return dz, ddz
  end
end

evalderiv(zeta::Vector{Complex128},map::ExteriorMap) = evalderiv(zeta,map,false)

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
            argr = argl + dist*(arg2-argl)
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
