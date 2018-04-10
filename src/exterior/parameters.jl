using NLsolve

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

    depfun! = Depfun(beta,nmlen,qdat)

    F0 = similar(y0)
    df = OnceDifferentiable(depfun!, y0, F0)
    sol = nlsolve(df,y0,show_trace = :false)

    zeta = zeros(Complex128,n)
    θ = zeros(n-1)
    y_to_zeta!(zeta,θ,sol.zero)

  end

  mid = zeta[1]*exp(0.5*im*angle(zeta[2]/zeta[1]))
  dequad = DQuad(beta,qdat)
  c = (w[2]-w[1])/(dequad([zeta[1]],[mid],[1],zeta)[1]-dequad([zeta[2]],[mid],[2],zeta)[1])

  return zeta, c

end

function evaluate_exterior(zeta::Vector{Complex128},w::Vector{Complex128},
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
  dequad = DQuad(beta,qdat)

  # initialize the mapped evaluation points
  z = zeros(Complex128,neval)

  # find the closest prevertices to each evaluation point and their
  #  corresponding distances
  dz = abs.(zeta[:,ones(Int,n)]-prev[:,ones(Int,neval)].')
  #dz = abs.(hcat([zeta for i=1:n]...)-vcat([transpose(prev) for i=1:neval]...))
  (dist,ind) = findmin(dz,2)
  sing = floor.(Int,(ind[:]-1)/neval)+1

  # find any prevertices in the evaluation list and set them equal to
  #  the corresponding vertices. The origin is also a singular point
  vertex = (dist[:] .< tol)
  z[vertex] = w[sing[vertex]]
  zerop = abs.(zeta) .< tol
  z[zerop] = Inf
  vertex = vertex .| zerop

  # the starting (closest) singularities for each evaluation point
  prevs = prev[sing]

  # set the initial values of the non-vertices
  z[.!vertex] = w[sing[.!vertex]]

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
    z[unf] = z[unf] + c*dequad(zetaold[unf],zetanew[unf],sing[unf],prev)

    # set new starting integration points for those that can be integrated
    #  further
    unf = dist .< 1
    zetaold[unf] = zetanew[unf]

    # only the first step can have a singularity
    sing .= 0

  end

  return z

end

function evalderiv_exterior(zeta::Vector{Complex128},beta::Vector{Float64},
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

function evalinv_exterior(z::Vector{Complex128},w::Vector{Complex128},
                  beta::Vector{Float64},prev::Vector{Complex128},
                  c::Complex128,qdat::Tuple{Array{Float64,2},Array{Float64,2}})

   n = length(w)
   zeta = zeros(Complex128,size(z))
   lenz = length(z)

   # Find z values close to vertices and set zeta to the corresponding
   # prevertices
   done = zeros(Bool,size(z))
   for j = 1:n
     idx = find(abs.(z-w[j]) .< 3*eps())
     zeta[idx] = prev[j]
     done[idx] = true
   end
   lenz -= sum(done)
   if lenz==0
     return zeta
   end

   # Now, for remaining z values, first try to integrate
   #  dζ/dt = (z - z(ζ₀))/z'(ζ) from t = 0 to t = 1,
   # with the initial condition ζ(0) = ζ₀.


   return zeta


#
# % ODE
# if ode
#   if isempty(zeta0)
#     % Pick a value z0 (not a singularity) and compute the map there.
#     [z0,w0] = scimapz0('de',z(~done),w,beta,prev,c,qdat);
#   else
#     w0 = demap(zeta0,w,beta,prev,c,qdat);
#     if length(zeta0)==1 & lenz > 1
#       zeta = zeta0(:,ones(lenz,1)).';
#       w0 = w0(:,ones(lenz,1)).';
#     end
#     w0 = w0(~done);
#     zeta0 = zeta0(~done);
#   end
#
#   % Use relaxed ODE tol if improving with Newton.
#   odetol = max(tol,1e-3*(newton));
#
#   % Rescale dependent coordinate
#   scale = (z(~done) - w0(:));
#
#   % Solve ODE
#   zeta0 = [real(zeta0);imag(zeta0)];
#   [t,y] = ode23('deimapfun',[0,0.5,1],zeta0,odeset('abstol',odetol),...
#       scale,prev,beta,c);
#   [m,leny] = size(y);
#   zeta(~done) = y(m,1:lenz)+sqrt(-1)*y(m,lenz+1:leny);
#   out = abs(zeta) > 1;
#   zeta(out) = sign(zeta(out));
# end
#
# % Newton iterations
# if newton
#   if ~ode
#     zetan = zeta0(:);
#     if length(zeta0)==1 & lenz > 1
#       zetan = zetan(:,ones(lenz,1));
#     end
#     zetan(done) = zeta(done);
#   else
#     zetan = zetap(:);
#   end
#
#   z = z(:);
#   k = 0;
#   while ~all(done) & k < maxiter
#     F = z(~done) - demap(zetan(~done),w,beta,prev,c,qdat);
#     m = length(F);
#     dF = c*(zetan(~done).').^(-2) .* exp(sum(beta(:,ones(m,1)) .* ...
#       log(1-(zetan(~done,ones(n,1)).')./prev(:,ones(m,1)))));
#     zetan(~done) = zetan(~done) + F(:)./dF(:);
#     done(~done) = (abs(F)< tol);
#     k = k+1;
#   end
#   if any(abs(F)> tol)
#     str = sprintf('Check solution; maximum residual = %.3g\n',max(abs(F)));
#     warning(str)
#   end
#   zeta(:) = zetan;
# end

end

function initial_guess(z::Vector{Complex128},w::Vector{Complex128},
                  beta::Vector{Float64},prev::Vector{Complex128},
                  c::Complex128,qdat::Tuple{Array{Float64,2},Array{Float64,2}})

  n = length(w)
  tol = 1000.0*10.0^(-size(qdat[1])[1])
  shape = copy(z)
  zeta0 = copy(z)
  z0 = copy(z)
  argz = angle.(prev);
  argz[argz.<=0] .+= 2π

  argw = cumsum([angle(w[2]-w[1]); -π*beta[2:n]])

  infty = isinf.(w)
  fwd = circshift(1:n,-1)
  anchor = zeros(w)
  anchor[.~infty] = w[.~infty]
  anchor[infty] = w[fwd[infty]]
  direcn = exp.(im*argw)
  direcn[infty] = -direcn[infty]
  len = abs.(w[fwd] - w)

  factor = 0.5
  done = zeros(Bool,length(z))
  M = length(z)
  iter = Int(0)

  A = zeros(Float64,2,2)

  zetabase = NaN*ones(Complex128,n)
  zbase = NaN*ones(Complex128,n)
  idx = []
  while M > 0
    for j = 1:n
      if j<n
        zetabase[j] = exp(im*(factor*argz[j] + (1-factor)*argz[j+1]))
      else
        zetabase[j] = exp(im*(factor*argz[n] + (1-factor)*(2π+argz[1])))
      end
    end
    zbase = evaluate_exterior(zetabase,w,beta,prev,c,qdat)
    proj = real.( (zbase-anchor) .* conj.(direcn) )
    zbase = anchor + proj.*direcn
    if isempty(idx)
      dist,idxtemp = findmin(abs.( z[.~done,ones(Int,n)].' - zbase[:,ones(Int,M)] ),1)
      for k = 1:M
          push!(idx,idxtemp[k]-(k-1)*n)
      end
    else
      idx[.~done] = idx[.~done].%n + 1
    end
    zeta0[.~done] = zetabase[idx[.~done]]
    z0[.~done] = zbase[idx[.~done]]

    for j = 1:n
      active = (idx.==j) .& (.~done)
      if any(active)

        done[active] = ones(Bool,sum(active))
        for k in [1:j-1;j+1:n]'
          A[:,1] = [real(direcn[k]);imag(direcn[k])]
          for p in find(active)
            dif = z0[p]-z[p]
              A[:,2] = [real(dif);imag(dif)]
              if cond(A) < eps()
                zx = real( (z[p]-anchor[k]) / direcn[k] )
                z0x = real( (z0[p]-anchor[k]) / direcn[k] )
                if (zx*z0x < 0) || ((zx-len[k])*(z0x-len[k]) < 0)
                  done[p] = false
                end
              else
                dif = z0[p]-anchor[k]
                  s = A\[real(dif);imag(dif)]
                  if s[1]>=0 && s[1]<=len[k]
                    if abs(s[2]-1) < tol
                      zeta0[p] = zetabase[k]
                      z0[p] = zbase[k]
                    elseif abs(s[2]) < tol
                      if real( conj(z[p]-z0[p])*im*direcn[k] ) > 0
                        done[p] = false
                      end
                    elseif s[2] > 0 && s[2] < 1
                      done[p] = false
                    end
                  end
                end # cond(A)
              end # for p
            end # for k
            M = sum(.~done)
            if M == 0
              break
            end
          end # if any active
          if iter > 2*n
            error("Can''t seem to choose starting points.  Supply them manually.")
          else
            iter += 1
          end
          factor = rand(1)[1]
        end # for j

  end  # while M

  return zeta0, z0

end

####

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

#######

struct DabsQuad{T,N,NQ}
  beta :: Vector{T}
  qdat :: Tuple{Array{T,2},Array{T,2}}
end

function DabsQuad(beta::Vector{T},tol::T) where T
  n = length(beta)
  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)
  DabsQuad{T,n,nqpts}(beta,qdat)
end

function DabsQuad(beta::Vector{T},
                  qdat::Tuple{Array{T,2},Array{T,2}}) where T
  n = length(beta)
  nqpts = size(qdat[1],1)
  DabsQuad{T,n,nqpts}(beta,qdat)
end

function (I::DabsQuad{T,N,NQ})(zeta1::Vector{Complex128},zeta2::Vector{Complex128},sing1::Vector{Int64},zeta::Vector{Complex128}) where {T,N,NQ}

   (qnode,qwght) = I.qdat
   argz = angle.(zeta)

   argz1 = angle.(zeta1)
   argz2 = angle.(zeta2)
   ang21 = angle.(zeta2./zeta1)

   bigargz = transpose(argz[:,ones(Int,NQ)])

   discont = (argz2-argz1).*ang21 .< 0
   argz2[discont] += 2π*sign.(ang21[discont])

   if isempty(sing1)
     sing1 = zeros(size(zeta1))
   end
   result = zeros(Float64,size(zeta1))

   nontriv = find(zeta1.!=zeta2)
    #tic()
   for k in nontriv
     z1k, z2k, arg1k, arg2k, sing1k =
          zeta1[k], zeta2[k], argz1[k], argz2[k], sing1[k]
     zetas = vcat(zeta[1:sing1k-1],zeta[sing1k+1:end])
     dist = min(1,2*minimum(abs.(zetas-z1k))/abs(z2k-z1k))
     argr = arg1k + dist*(arg2k-arg1k)
     ind = ((sing1k+N) % (N+1)) + 1
     nd = 0.5*((argr-arg1k)*qnode[:,ind] + argr + arg1k)
     wt = 0.5*abs(argr-arg1k)*qwght[:,ind]
     θ = (nd[:,ones(Int,N)]-bigargz.+2π).%(2π)
     θ[θ.>π] = 2π-θ[θ.>π]
     terms = 2sin.(0.5θ)
     if !any(terms==0.0)
        if sing1k > 0
            terms[:,sing1k] ./= abs.(nd-arg1k)
            wt .*= (0.5*abs.(argr-arg1k)).^I.beta[sing1k]
        end
        result[k] = At_mul_B(exp.(log.(terms)*I.beta),wt)
        while dist < 1.0
            argl = argr
            zetal = exp(im*argl)
            dist = min(1,2*minimum(abs.(zeta-zetal)/abs(zetal-z2k)))
            argr = argl + dist*(arg2k-argl)
            nd = 0.5*((argr-argl)*qnode[:,N+1] + argr + argl)
            wt = 0.5*abs(argr-argl)*qwght[:,N+1]
            #θ = hcat([(nd - argz[i] + 2π).%(2π) for i = 1:N]...)
            θ = (nd[:,ones(Int,N)]-bigargz.+2π).%(2π)
            θ[θ.>π] = 2π-θ[θ.>π]
            terms = 2sin.(0.5θ)
            result[k] += At_mul_B(exp.(log.(terms)*I.beta),wt)
        end
     end
   end
    #toc()
    #nothing
   return result
end

struct DQuad{T,N,NQ}
  beta :: Vector{T}
  qdat :: Tuple{Array{T,2},Array{T,2}}
end

function DQuad(beta::Vector{T},tol::T) where T
  n = length(beta)
  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(beta,nqpts)
  DQuad{T,n,nqpts}(beta,qdat)
end

function DQuad(beta::Vector{T},
                  qdat::Tuple{Array{T,2},Array{T,2}}) where T
  n = length(beta)
  nqpts = size(qdat[1],1)
  DQuad{T,n,nqpts}(beta,qdat)
end


function (I::DQuad{T,N,NQ})(zeta1::Vector{Complex128},zeta2::Vector{Complex128},
          sing1::Vector{Int64},zeta::Vector{Complex128}) where {T,N,NQ}

   (qnode,qwght) = I.qdat

   beta = [I.beta;-2]

   bigzeta = transpose(zeta[:,ones(Int,NQ)])

   if isempty(sing1)
     sing1 = zeros(Int,size(zeta1))
   end
   result = zeros(Complex128,size(zeta1))

   nontriv = find(zeta1.!=zeta2)
    #tic()
   for k in nontriv
     z1k, z2k, sing1k = zeta1[k], zeta2[k], sing1[k]
     zetas = vcat(zeta[1:sing1k-1],zeta[sing1k+1:end])
     dist = min(1,2*minimum(abs.(zetas-z1k))/abs(z2k-z1k))
     zetar = z1k + dist*(z2k-z1k)
     ind = sing1k + (N+1)*(sing1k==0)

     nd = 0.5*((zetar-z1k)*qnode[:,ind] + zetar + z1k)
     wt = 0.5*(zetar-z1k)*qwght[:,ind]
     terms = 1 .- nd[:,ones(Int,N)]./bigzeta
     if !any(terms==0.0)
       terms = hcat(terms,nd)
        if sing1k > 0
            terms[:,sing1k] ./= abs.(terms[:,sing1k])
            wt .*= (0.5*abs.(zetar-z1k)).^beta[sing1k]
        end
        result[k] = At_mul_B(exp.(log.(terms)*beta),wt)

        while dist < 1.0
            zetal = zetar
            dist = min(1,2*minimum(abs.(zeta-zetal)/abs(zetal-z2k)))
            zetar = zetal + dist*(z2k-zetal)
            nd = 0.5*((zetar-zetal)*qnode[:,N+1] + zetar + zetal)
            wt = 0.5*(zetar-zetal)*qwght[:,N+1]
            terms = 1 .- nd[:,ones(Int,N)]./bigzeta
            terms = hcat(terms,nd)
            result[k] += At_mul_B(exp.(log.(terms)*beta),wt)
        end
      end
    end
   return result
end


struct Depfun{T,N,NQ}
    zeta :: Vector{Complex{T}}
    beta :: Vector{T}
    nmlen :: Vector{T}
    qdat :: Tuple{Array{T,2},Array{T,2}}
    θ    :: Vector{T}
    mid  :: Vector{Complex{T}}
    ints :: Vector{T}
    dabsquad :: DabsQuad{T,N,NQ}
end

function Depfun(beta::Vector{T},nmlen::Vector{T},qdat:: Tuple{Array{T,2},Array{T,2}}) where T
    # should compute nmlen in here
    n = length(beta)
    nqpts = size(qdat[1],1)
    zeta = zeros(Complex128,n)
    θ = zeros(n-1)
    mid = zeros(Complex128,n-2)
    ints = zeros(Complex128,n-2)
    dabsquad = DabsQuad(beta,qdat)
    Depfun{T,n,nqpts}(zeta,beta,nmlen,qdat,θ,mid,ints,dabsquad)
end

function (R::Depfun{T,N,NQ})(F,y) where {T,N,NQ}


  y_to_zeta!(R.zeta,R.θ,y)

  @. R.mid = exp(im*0.5*(R.θ[1:N-2]+R.θ[2:N-1]))

  #tic()
  R.ints .= R.dabsquad(R.zeta[1:N-2],R.mid,collect(1:N-2),R.zeta)
  R.ints .+= R.dabsquad(R.zeta[2:N-1],R.mid,collect(2:N-1),R.zeta)

  #toc()

  if N > 3
    @. F[1:N-3] = abs(R.ints[2:N-2])/abs(R.ints[1]) - R.nmlen
  end

  res = -sum(R.beta./R.zeta)/R.ints[1]
  @. F[N-2] = real(res)
  @. F[N-1] = imag(res)
end


function y_to_zeta!(zeta::Vector{Complex{T}},
                    θ::Vector{T},y::Vector{T}) where T

  cs = cumsum(cumprod([1;exp.(-y)]))
  n = length(cs)
  @. θ = 2π*cs[1:n-1]/cs[n]
  zeta[n] = 1.0
  @. zeta[1:n-1] = exp(im*θ)

end
