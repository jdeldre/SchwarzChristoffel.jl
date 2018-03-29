module Exterior

using NLsolve
using ..Properties
using ..Polygons
using ..Integration

export Extermap,dabsquad

struct Extermap <: Map

  z    :: Vector{Complex128}
  beta :: Vector{Float64}
  qdat :: Tuple{Array{Float64,2},Array{Float64,2}}
  tol  :: Float64

  zeta :: Vector{Complex128}
  C    :: Complex128

end

function Extermap(p::Polygon)

  n = length(p.vert)

  zeta0 = Complex128[]
  tol = 1e-8

  w = flipdim(vertex(p),1)
  beta = 1.-flipdim(interiorangle(p),1)

  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)

  zeta, c = deparam(w,beta,zeta0)

  Extermap(p.vert,p.angle,qdat,tol,zeta,c)
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

function DabsQuad(zeta::Vector{Complex128},beta::Vector{Float64},qdat::Tuple{Array{Float64,2},Array{Float64,2}})
  n = length(zeta)
  nqpts = size(qdat[1],1)
  DabsQuad{n}(zeta,beta,nqpts,qdat)
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
   result = 0
   if zeta1 != zeta2
     dist = min(1,2*minimum(abs.([I.zeta[1:sing1-1];I.zeta[sing1+1:n]]-zeta1))/abs(zeta2-zeta1))
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
        result = dot(exp.(log.(terms)*I.beta),wt)
        while dist < 1.0
            argl = argr
            zl = exp(im*argl)
            dist = min(1,2*minimum(abs.(I.zeta-zl)/abs(zl-zeta2)))
            argr = argl + dist*(arg2-argl)
            nd = 0.5*((argr-argl)*qnode[:,n+1] + argr + argl)
            wt = 0.5*abs(argr-argl)*qwght[:,n+1]
            θ = hcat([(nd - argz[i] + 2π).%(2π) for i = 1:n]...)
            θ[θ.>π] = 2π-θ[θ.>π]
            result += dot(exp.(log.(2sin.*(0.5θ))*I.beta),wt)
        end
    end
   end
   return result
end


function deparam(w::Vector{Complex128},beta::Vector{Float64},
                 zeta0::Vector{Complex128},qdat)
  # w clockwise
  # beta turning angles

  n = length(w)
  if n == 2
    zeta = [-1,1]
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
    sol = nlsolve(depfun!,y0,show_trace = :true)

    zeta, θ = y_to_zeta(sol.zero)

  end

  mid = zeta[1]*exp(0.5*im*angle(zeta[2]/zeta[1]))
  dabsquad = DabsQuad(zeta,beta,qdat)
  c = 0.0
  #c = (w[2]-w[1])/(dequad(zeta[1],mid,1)-dequad(zeta[2],mid,2))

  return zeta, c

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
