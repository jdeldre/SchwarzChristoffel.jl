#=
Many of the base level routines here are adapted from the SC Toolbox by Toby Driscoll
and are inspired from the work of L. Trefethen, "Numerical Computation of the
Schwarz-Christoffel Transformation", STAN-CS-79-710, 1979.
=#


#= functions for the Schwarz-Christoffel exterior map from unit disk =#
function param(w::Vector{ComplexF64},β::Vector{Float64},
                 ζ0::Vector{ComplexF64},
                 qdat::Tuple{Array{Float64,2},Array{Float64,2}})
#=
Solve for the parameters of the exterior Schwarz-Christoffel mapping: the
prevertices `ζ` and the constant factor `c`. The routine requires the
vertices `w` in clockwise order, the exterior turning angles `β`,
a list of guesses for the prevertices `ζ0` (which can be blank), and
the Gauss-Jacobi quadrature data in `qdat`
=#

  n = length(w)
  if n == 2
    ζ = ComplexF64[-1,1]
  else
    len = abs.(diff(circshift(w,1)))
    nmlen = abs.(len[3:n-1]/len[2])
    if isempty(ζ0)
      y0 = zeros(n-1)
    else
      ζ0 = ζ0/ζ0[n]
      θ = angle.(ζ0)
      θ[θ.<=0] = θ[θ.<=0] .+ 2π
      dt = diff([0;θ[1:n-1];2π])
      @. y0 = log(dt[1:n-1]/dt[2:n])
    end

    depfun! = Depfun(β,nmlen,qdat)

    F0 = similar(y0)
    df = OnceDifferentiable(depfun!, y0, F0)
    sol = nlsolve(df,y0,show_trace = :false)

    ζ = zeros(ComplexF64,n)
    θ = zeros(n-1)
    y_to_ζ!(ζ,θ,sol.zero)

  end

  mid = ζ[1]*exp(0.5*im*angle(ζ[2]/ζ[1]))
  dequad = DQuad(β,qdat)
  c = (w[2]-w[1])/(dequad([ζ[1]],[mid],[1],ζ)[1]-dequad([ζ[2]],[mid],[2],ζ)[1])

  return ζ, c

end

function evaluate_exterior(ζ::Vector{ComplexF64},w::Vector{ComplexF64},
                  β::Vector{Float64},prev::Vector{ComplexF64},
                  c::ComplexF64,qdat::Tuple{Array{Float64,2},Array{Float64,2}})
  #=
  Evaluates the exterior Schwarz-Christoffel mapping at `ζ`, which is
  presumed to be inside the unit circle. The vector `w` are the vertices (in
  clockwise order), `β` are the exterior turning angles (also in clockwise order),
  `prev` are the prevertices on the unit circle, and `c` is the constant factor
  of the transformation. The tuple `qdat` contains the Gauss-Jacobi nodes and
  weights.
  =#

  if isempty(ζ)
    nothing
  end

  n = length(w)
  neval = length(ζ)
  tol = 10.0^(-size(qdat[1],1))

  # set up the integrator
  dequad = DQuad(β,qdat)

  # initialize the mapped evaluation points
  z = zeros(ComplexF64,neval)

  # find the closest prevertices to each evaluation point and their
  #  corresponding distances
  dz = abs.(repeat(ζ, 1, n) - repeat(transpose(prev), neval))
  #dz = abs.(hcat([ζ for i=1:n]...)-vcat([transpose(prev) for i=1:neval]...))
  (dist,ind) = findmin(dz, dims = 2)
  #sing = floor.(Int,(ind[:]-1)/neval)+1
  sing = [I[2] for I in ind]

  # find any prevertices in the evaluation list and set them equal to
  #  the corresponding vertices. The origin is also a singular point
  vertex = (dist[:] .< tol)
  z[vertex] = w[sing[vertex]]
  zerop = abs.(ζ) .< tol
  z[zerop] .= Inf
  vertex = vertex .| zerop

  # the starting (closest) singularities for each evaluation point
  prevs = prev[sing]

  # set the initial values of the non-vertices
  z[.!vertex] = w[sing[.!vertex]]

  # distance to singularity at the origin
  absζ = abs.(ζ)

  # unfinished cases
  unf = .!vertex

  # set the integration starting points
  ζold = copy(prevs)
  ζnew = copy(ζold)
  dist = ones(neval)
  while any(unf)
    # distance to integrate still
    dist[unf] = min.(1,2*absζ[unf]./abs.(ζ[unf]-ζold[unf]))

    # new integration end point
    ζnew[unf] = ζold[unf] + dist[unf].*(ζ[unf]-ζold[unf])

    # integrate
    z[unf] = z[unf] + c*dequad(ζold[unf],ζnew[unf],sing[unf],prev)

    # set new starting integration points for those that can be integrated
    #  further
    unf = dist .< 1
    ζold[unf] = ζnew[unf]

    # only the first step can have a singularity
    sing .= 0

  end

  return z

end

function evalderiv_exterior(ζ::Vector{ComplexF64},β::Vector{Float64},
                  prev::Vector{ComplexF64},c::ComplexF64)

#=
Evaluates the first and second derivative of the exterior Schwarz-Christoffel
mapping at `ζ`, which is
presumed to be inside the unit circle. The vector `β` are the exterior turning
angles (also in clockwise order), `prev` are the prevertices on the unit circle, and
`c` is the constant factor of the transformation.
=#
    n = length(prev)
    neval = length(ζ)
    β = [β;-2]
    terms = [hcat([1 .- ζ/prev[i] for i = 1:n]...) ζ]
    dz = c*exp.(log.(terms)*β)
    terms2 = [[-β[i]/prev[i] for i = 1:n];β[n+1]]
    ddz = dz.*((1.0./terms)*terms2)
    return dz, ddz

end

function evalinv_exterior(z::Vector{ComplexF64},w::Vector{ComplexF64},
                  β::Vector{Float64},prev::Vector{ComplexF64},
                  c::ComplexF64,qdat::Tuple{Array{Float64,2},Array{Float64,2}})

#=
Evaluates the inverse of the exterior Schwarz-Christoffel mapping, using a combination
of integration and Newton iteration, using techniques from Trefethen (1979).
=#

   n = length(w)
   ζ = zeros(ComplexF64,size(z))
   lenz = length(z)
   ζ0 = []
   maxiter = 10
   tol = 10.0^(-size(qdat[1],1))

   # Find z values close to vertices and set ζ to the corresponding
   # prevertices
   done = zeros(Bool,size(z))
   for j = 1:n
     idx = findall(abs.(z .- w[j]) .< 3*eps())
     ζ[idx] .= prev[j]
     done[idx] .= true
   end
   lenz -= sum(done)
   if lenz==0
     return ζ
   end

   # Now, for remaining z values, first try to integrate
   #  dζ/dt = (z - z(ζ₀))/z'(ζ) from t = 0 to t = 1,
   # with the initial condition ζ(0) = ζ₀.
   if isempty(ζ0)
     ζ0,z0 = initial_guess(z,w,β,prev,c,qdat)
   else
     z0 = evaluate_exterior(ζ0,w,β,prev,c,qdat)
     if length(ζ0)==1 && lenz > 1
       ζ0 = repeat(transpose(ζ0), lenz)
       z0 = repeat(transpose(z0), lenz)
     end
     z0 = z0[.~done]
     ζ0 = ζ0[.~done]
   end
   odetol = max(tol,1e-3)
   scale = z[.~done] - z0

   ζ0 = [real(ζ0);imag(ζ0)]

   f(ζ,p,t) = invfunc(ζ,scale,β,prev,c)
   tspan = (0.0,1.0)
   prob = ODEProblem(f,ζ0,tspan)
   sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
   lenu = length(ζ0)
   ζ[.~done] = sol.u[end][1:lenz]+im*sol.u[end][lenz+1:lenu];
   out = abs.(ζ) .> 1
   ζ[out] = sign.(ζ[out])

   # Now use Newton iterations to improve the solution
   ζn = ζ
   k = 0
   while ~all(done) && k < maxiter
     F = z[.~done] - evaluate_exterior(ζn[.~done],w,β,prev,c,qdat)
     M = length(F)
     dF, ddz = evalderiv_exterior(ζn[.~done],β,prev,c)
     ζn[.~done] = ζn[.~done] + F./dF

     done[.~done] = abs.(F).< tol
     k += 1
   end
   F = z[.~done] - evaluate_exterior(ζn[.~done],w,β,prev,c,qdat)
   if any(abs.(F).> tol)
     error("Check solution")
   end
   ζ = ζn

end


function invfunc(u,scale,β::Vector{Float64},prev::Vector{ComplexF64},c::ComplexF64)
    lenu = length(u)
    lenzp = Int(lenu/2)
    ζ = u[1:lenzp]+im*u[lenzp+1:lenu]

    dz, ddz = evalderiv_exterior(ζ,β,prev,c)
    f = scale./dz
    zdot = [real(f);imag(f)]
end

function initial_guess(z::Vector{ComplexF64},w::Vector{ComplexF64},
                  β::Vector{Float64},prev::Vector{ComplexF64},
                  c::ComplexF64,qdat::Tuple{Array{Float64,2},Array{Float64,2}})

  n = length(w)
  tol = 1000.0*10.0^(-size(qdat[1])[1])
  shape = copy(z)
  ζ0 = copy(z)
  z0 = copy(z)
  argz = angle.(prev);
  argz[argz.<=0] .+= 2π

  argw = cumsum([angle(w[2]-w[1]); -π*β[2:n]])

  infty = isinf.(w)
  fwd = circshift(1:n,-1)
  anchor = zero(w)
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

  ζbase = NaN*ones(ComplexF64,n)
  zbase = NaN*ones(ComplexF64,n)
  idx = []
  while M > 0
    for j = 1:n
      if j<n
        ζbase[j] = exp(im*(factor*argz[j] + (1-factor)*argz[j+1]))
      else
        ζbase[j] = exp(im*(factor*argz[n] + (1-factor)*(2π+argz[1])))
      end
    end
    zbase = evaluate_exterior(ζbase,w,β,prev,c,qdat)
    proj = real.( (zbase-anchor) .* conj.(direcn) )
    zbase = anchor + proj.*direcn
    if isempty(idx)
      dist,idxtemp = findmin(abs.( repeat(transpose(z[.~done]), n) - repeat(zbase, 1, M)), dims = 1)
      for k = 1:M
          push!(idx,idxtemp[k][1])
      end
    else
      idx[.~done] = idx[.~done].%n .+ 1
    end
    ζ0[.~done] = ζbase[idx[.~done]]
    z0[.~done] = zbase[idx[.~done]]

    for j = 1:n
      active = (idx.==j) .& (.~done)
      if any(active)

        done[active] = ones(Bool,sum(active))
        for k in [1:j-1;j+1:n]'
          A[:,1] = [real(direcn[k]);imag(direcn[k])]
          for p in findall(active)
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
                      ζ0[p] = ζbase[k]
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

  return ζ0, z0

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

  # Find the coefficient for 1/ζ^(pow-1) in the fhat expansion
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
  Pplus1 = sortslices(Pplus1, dims=1)
  dP = sum(abs.([transpose(Pplus1[1,:]);diff(Pplus1,1)]),2)
  return Pplus1[dP[:].!=0,:]

end

#######

struct DabsQuad{T,N,NQ}
  β :: Vector{T}
  qdat :: Tuple{Array{T,2},Array{T,2}}
end

function DabsQuad(β::Vector{T},tol::T) where T
  n = length(β)
  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(β,nqpts)
  DabsQuad{T,n,nqpts}(β,qdat)
end

function DabsQuad(β::Vector{T},
                  qdat::Tuple{Array{T,2},Array{T,2}}) where T
  n = length(β)
  nqpts = size(qdat[1],1)
  DabsQuad{T,n,nqpts}(β,qdat)
end

function (I::DabsQuad{T,N,NQ})(ζ1::Vector{ComplexF64},ζ2::Vector{ComplexF64},
                        sing1::Vector{Int64},ζ::Vector{ComplexF64}) where {T,N,NQ}
#=
This integrates |z'(λ)| from `ζ1` to `ζ2` on the unit circle, using
Gauss-Jacobi quadrature, where
`sing1` contains the index of the prevertex in `ζ` that `ζ1` corresponds
to, or 0 if `ζ1` is not a prevertex. Note that `ζ2` cannot be a prevertex.
=#

   (qnode,qwght) = I.qdat

   argz = angle.(ζ)

   argζ1 = angle.(ζ1)
   argζ2 = angle.(ζ2)
   ang21 = angle.(ζ2./ζ1)

   bigargz = transpose(repeat(argz, 1, NQ))

   discont = (argζ2-argζ1).*ang21 .< 0
   argζ2[discont] += 2π*sign.(ang21[discont])

   if isempty(sing1)
     sing1 = zeros(size(ζ1))
   end
   result = zeros(Float64,size(ζ1))

   nontriv = findall(ζ1.!=ζ2)
    #tic()
   for k in nontriv
     ζ1k, ζ2k, arg1k, arg2k, sing1k =
          ζ1[k], ζ2[k], argζ1[k], argζ2[k], sing1[k]
     ζs = vcat(ζ[1:sing1k-1],ζ[sing1k+1:end])
     dist = min(1,2*minimum(abs.(ζs .- ζ1k))/abs(ζ2k-ζ1k))
     argr = arg1k + dist*(arg2k-arg1k)
     ind = ((sing1k+N) % (N+1)) + 1
     ζnd = 0.5*((argr-arg1k)*qnode[:,ind] .+ argr .+ arg1k)
     wt = 0.5*abs(argr-arg1k)*qwght[:,ind]
     θ = (repeat(ζnd, 1, N) .- bigargz .+ 2π).%(2π)
     θ[θ.>π] = 2π .- θ[θ.>π]
     terms = 2sin.(0.5θ)
     if !any(terms==0.0)
        if sing1k > 0
            terms[:,sing1k] ./= abs.(ζnd .- arg1k)
            wt .*= (0.5*abs.(argr-arg1k)).^I.β[sing1k]
        end
        result[k] = transpose(exp.(log.(terms)*I.β))*wt
        while dist < 1.0
            argl = argr
            ζl = exp(im*argl)
            dist = min(1,2*minimum(abs.(ζ .- ζl)/abs(ζl-ζ2k)))
            argr = argl + dist*(arg2k-argl)
            ζnd = 0.5*((argr-argl)*qnode[:,N+1] .+ argr .+ argl)
            wt = 0.5*abs(argr-argl)*qwght[:,N+1]
            #θ = hcat([(ζnd - argz[i] + 2π).%(2π) for i = 1:N]...)
            θ = (repeat(ζnd, 1, N) .- bigargz .+ 2π).%(2π)
            θ[θ.>π] = 2π .- θ[θ.>π]
            terms = 2sin.(0.5θ)
            result[k] += transpose(exp.(log.(terms)*I.β))*wt
        end
     end
   end
    #toc()
    #nothing
   return result
end

function (I::DabsQuad{T,N,NQ})(ζ1::Vector{ComplexF64},ζ2::Vector{ComplexF64},
                        sing1::Vector{Int64},ζ::Vector{ComplexF64},pow::Int) where {T,N,NQ}

#=
This integrates λ⁻ᵏz'(λ) (where k is `pow`) from `ζ1` to `ζ2` on the unit circle,
using Gauss-Jacobi quadrature, where
`sing1` contains the index of the prevertex in `ζ` that `ζ1` corresponds
to, or 0 if `ζ1` is not a prevertex. Note that `ζ2` cannot be a prevertex.
=#

   (qnode,qwght) = I.qdat

   argz = angle.(ζ)

   argζ1 = angle.(ζ1)
   argζ2 = angle.(ζ2)
   ang21 = angle.(ζ2./ζ1)

   bigargz = transpose(repeat(argz, 1, NQ))

   discont = (argζ2-argζ1).*ang21 .< 0
   argζ2[discont] += 2π*sign.(ang21[discont])

   if isempty(sing1)
     sing1 = zeros(size(ζ1))
   end
   result = zeros(ComplexF64,size(ζ1))

   nontriv = findall(ζ1.!=ζ2)
    #tic()
   for k in nontriv
     ζ1k, ζ2k, arg1k, arg2k, sing1k =
          ζ1[k], ζ2[k], argζ1[k], argζ2[k], sing1[k]
     ζs = vcat(ζ[1:sing1k-1],ζ[sing1k+1:end])
     dist = min(1,2*minimum(abs.(ζs .- ζ1k))/abs(ζ2k-ζ1k))
     argr = arg1k + dist*(arg2k-arg1k)
     ind = ((sing1k+N) % (N+1)) + 1
     ζnd = 0.5*((argr-arg1k)*qnode[:,ind] .+ argr .+ arg1k)
     wt = 0.5*abs(argr-arg1k)*qwght[:,ind]
     θ = (repeat(ζnd, 1, N) .- bigargz .+ 2π).%(2π)
     #θ[θ.>π] = 2π-θ[θ.>π]
     #terms = 2sin.(0.5θ)
     terms = 1 .- exp.(im*θ)
     if !any(terms==0.0)
        #terms = hcat(terms,exp.(im*ζnd))
        if sing1k > 0
            terms[:,sing1k] ./= abs.(ζnd .- arg1k)
            wt .*= (0.5*abs.(argr-arg1k)).^I.β[sing1k]
        end
        result[k] = transpose(exp.(log.(terms)*I.β+im*(pow-1)*ζnd))*wt
        while dist < 1.0
            argl = argr
            ζl = exp(im*argl)
            dist = min(1,2*minimum(abs.(ζ .- ζl)/abs(ζl-ζ2k)))
            argr = argl + dist*(arg2k-argl)
            ζnd = 0.5*((argr-argl)*qnode[:,N+1] .+ argr .+ argl)
            wt = 0.5*abs(argr-argl)*qwght[:,N+1]
            #θ = hcat([(ζnd - argz[i] + 2π).%(2π) for i = 1:N]...)
            θ = (repeat(ζnd, 1, N) .- bigargz .+ 2π).%(2π)
            #θ[θ.>π] = 2π-θ[θ.>π]
            #terms = 2sin.(0.5θ)
            terms = 1 .- exp.(im*θ)
            #terms = hcat(terms,exp.(im*ζnd))
            result[k] += transpose(exp.(log.(terms)*I.β+im*(pow-1)*ζnd))*wt
        end
     end
   end
    #toc()
    #nothing
   return result
end

struct DQuad{T,N,NQ}
  β :: Vector{T}
  qdat :: Tuple{Array{T,2},Array{T,2}}
end

function DQuad(β::Vector{T},tol::T) where T
  n = length(β)
  nqpts = max(ceil(Int,-log10(tol)),2)
  qdat = qdata(β,nqpts)
  DQuad{T,n,nqpts}(β,qdat)
end

function DQuad(β::Vector{T},
                  qdat::Tuple{Array{T,2},Array{T,2}}) where T
  n = length(β)
  nqpts = size(qdat[1],1)
  DQuad{T,n,nqpts}(β,qdat)
end


function (I::DQuad{T,N,NQ})(ζ1::Vector{ComplexF64},ζ2::Vector{ComplexF64},
          sing1::Vector{Int64},ζ::Vector{ComplexF64};pow::Int=0) where {T,N,NQ}
#=
This integrates z'(λ) from `ζ1` to `ζ2` on a straight path in the circle
plane, where `sing1` contains the index of the prevertex in `ζ` that `ζ1` corresponds
to, or 0 if `ζ1` is not a prevertex. Note that `ζ2` cannot be a prevertex.
If the optional argument `pow` is included, then it integrates λ⁻ᵏz'(λ), where
k = `pow`.
=#

   (qnode,qwght) = I.qdat

   β = [I.β;-2+pow]

   bigζ = transpose(repeat(ζ, 1, NQ))

   if isempty(sing1)
     sing1 = zeros(Int,size(ζ1))
   end
   result = zeros(ComplexF64,size(ζ1))

   nontriv = findall(ζ1.!=ζ2)
    #tic()
   for k in nontriv
     ζ1k, ζ2k, sing1k = ζ1[k], ζ2[k], sing1[k]
     ζs = vcat(ζ[1:sing1k-1],ζ[sing1k+1:end])
     dist = min(1,2*minimum(abs.(ζs .- ζ1k))/abs(ζ2k-ζ1k))
     ζr = ζ1k + dist*(ζ2k-ζ1k)
     # Choose which type of Gauss-Jacobi weights to use based on whether
     # ζ1k is a prevertex.
     ind = sing1k + (N+1)*(sing1k==0)

     ζnd = 0.5*((ζr-ζ1k)*qnode[:,ind] .+ ζr .+ ζ1k)
     wt = 0.5*(ζr-ζ1k)*qwght[:,ind]
     terms = 1 .- repeat(ζnd, 1, N)./bigζ
     if !any(terms==0.0)
       terms = hcat(terms,ζnd)
        # If ζ1k is a prevertex, adjust the integrand so that it is properly
        # set up for the Gauss-Jacobi integration
        if sing1k > 0
            terms[:,sing1k] ./= abs.(terms[:,sing1k])
            wt .*= (0.5*abs.(ζr-ζ1k)).^β[sing1k]
        end
        result[k] = transpose(exp.(log.(terms)*β))*wt

        while dist < 1.0
            ζl = ζr
            dist = min(1,2*minimum(abs.(ζ .- ζl)/abs(ζl-ζ2k)))
            ζr = ζl + dist*(ζ2k-ζl)
            ζnd = 0.5*((ζr-ζl)*qnode[:,N+1] .+ ζr .+ ζl)
            wt = 0.5*(ζr-ζl)*qwght[:,N+1]
            terms = 1 .- repeat(ζnd, 1, N)./bigζ
            terms = hcat(terms,ζnd)
            result[k] += transpose(exp.(log.(terms)*β))*wt
        end
      end
    end
   return result
end


struct Depfun{T,N,NQ}
    ζ :: Vector{Complex{T}}
    β :: Vector{T}
    nmlen :: Vector{T}
    qdat :: Tuple{Array{T,2},Array{T,2}}
    θ    :: Vector{T}
    mid  :: Vector{Complex{T}}
    ints :: Vector{T}
    dabsquad :: DabsQuad{T,N,NQ}
end

function Depfun(β::Vector{T},nmlen::Vector{T},qdat:: Tuple{Array{T,2},Array{T,2}}) where T
    # should compute nmlen in here
    n = length(β)
    nqpts = size(qdat[1],1)
    ζ = zeros(ComplexF64,n)
    θ = zeros(n-1)
    mid = zeros(ComplexF64,n-2)
    ints = zeros(ComplexF64,n-2)
    dabsquad = DabsQuad(β,qdat)
    Depfun{T,n,nqpts}(ζ,β,nmlen,qdat,θ,mid,ints,dabsquad)
end

function (R::Depfun{T,N,NQ})(F,y) where {T,N,NQ}


  y_to_ζ!(R.ζ,R.θ,y)

  @. R.mid = exp(im*0.5*(R.θ[1:N-2]+R.θ[2:N-1]))

  #tic()
  R.ints .= R.dabsquad(R.ζ[1:N-2],R.mid,collect(1:N-2),R.ζ)
  R.ints .+= R.dabsquad(R.ζ[2:N-1],R.mid,collect(2:N-1),R.ζ)

  #toc()

  if N > 3
    @. F[1:N-3] = abs(R.ints[2:N-2])/abs(R.ints[1]) - R.nmlen
  end

  res = -sum(R.β./R.ζ)/R.ints[1]
  F[N-2] = real(res)
  F[N-1] = imag(res)
end


function y_to_ζ!(ζ::Vector{Complex{T}},
                    θ::Vector{T},y::Vector{T}) where T

  cs = cumsum(cumprod([1;exp.(-y)]))
  n = length(cs)
  @. θ = 2π*cs[1:n-1]/cs[n]
  ζ[n] = 1.0
  @. ζ[1:n-1] = exp(im*θ)

end
