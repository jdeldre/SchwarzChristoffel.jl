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

end

function Extermap(p::Polygon)

  n = length(p.vert)

  zetai = Complex128[]
  tol = 1e-8

  w = flipdim(vertex(p),1)
  beta = 1.-flipdim(interiorangle(p),1)

  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)


  Extermap(p.vert,p.angle,qdat,tol)
end

function dabsquad(z1::Complex128,z2::Complex128,
                  sing1::Int64,z::Vector{Complex128},
                  beta::Vector{Float64},qdat)

   (qnode,qwght) = qdat
   nqpts = size(qnode,1)
   n = length(z)
   argz = angle.(z)

   argz1 = angle(z1)
   argz2 = angle(z2)
   ang21 = angle(z2/z1)

   if (argz2-argz1)*ang21 < 0
      argz2 += 2π*sign(ang21)
   end
   if isempty(sing1)
     sing1 = 0
   end
   I = 0
   if z1 != z2
     dist = min(1,2*minimum(abs.([z[1:sing1-1];z[sing1+1:n]]-z1))/abs(z2-z1))
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
            wt .*= (0.5*abs.(argr-argz1)).^beta[sing1]
        end
        I = dot(exp.(log.(terms)*beta),wt)
        while dist < 1.0
            argl = argr
            zl = exp(im*argl)
            dist = min(1,2*minimum(abs.(z-zl)/abs(zl-z2)))
            argr = argl + dist*(arg2-argl)
            nd = 0.5*((argr-argl)*qnode[:,n+1] + argr + argl)
            wt = 0.5*abs(argr-argl)*qwght[:,n+1]
            θ = hcat([(nd - argz[i] + 2π).%(2π) for i = 1:n]...)
            θ[θ.>π] = 2π-θ[θ.>π]
            I += dot(exp.(log.(2sin.*(0.5θ))*beta),wt)
        end
    end
   end
   return I
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

    F = zeros(Float64,n-1)
    depfun!(F,y) = depfunfull!(F,y,beta,nmlen,qdat)

    y = nlsolve(depfun!,y0,autodiff = :forward)
  end


end

function depfunfull!(F,y,beta,nmlen,qdat)

  cs = cumsum(cumprod([1;exp.(-y)]))
  θ = 2π*cs[1:n-1]./cs[n]
  zeta = ones(Complex128,n)
  @. zeta[1:n-1] = exp(im*θ)
  mid = zeros(Complex128,n-2)
  @. mid = exp(im*0.5*(θ[1:n-2]+θ[2:n-1]))

  dabsquad(z1,z2,sing1) = dabsquad(z1,z2,sing1,zeta,beta,qdat)

  ints = dabsquad.(zeta[1:n-2],mid,collect(1:n-2))
       + dabsquad.(zeta[2:n-1],mid,collect(2:n-1))

  if n > 3
    Ffill = abs.(ints[2:n-2])/abs(ints[1]) - nmlen
  end

  res = -sum(beta./zeta)/ints[1]
  F = [Ffill;real.(res);imag.(res)]
end

end
