#= functions for the power series exterior map from unit disk =#

struct PowerSeries{N,T}
  ccoeff :: Vector{T}
  dcoeff :: Vector{T}
end

function PowerSeries(ccoeff::Vector{T},dcoeff::Vector{T}) where T
  ncoeff = length(ccoeff)-2
  PowerSeries{ncoeff,T}(ccoeff,dcoeff)
end


function Base.show(io::IO,ps::PowerSeries{N,T}) where {N,T}
  println(io, "multipole coefficients:")
  println(io, "  c₁ = $(ps.ccoeff[1]), ")
  println(io, "  c₀ = $(ps.ccoeff[2]), ")
  print(io,   "  c₋ᵢ = ")
  for i = 1:N
    print(io,"$(ps.ccoeff[2+i]), ")
  end
  println(io, "i = 1:$(N)")
end

function (ps::PowerSeries)(ζ::T) where T<:Number
  ζⁿ = ζ
  z = zero(ζ)
  for c in ps.ccoeff
    z += c*ζⁿ
    ζⁿ /= ζ
  end
  return z
end

(ps::PowerSeries)(ζ::AbstractArray{T}) where T<:Number = ps.(ζ)

struct PowerSeriesDerivatives
  ps :: PowerSeries
end

function (dps::PowerSeriesDerivatives)(ζ::T) where T<:Number
  C = dps.ps.ccoeff
  dz = C[1]
  ζⁿ = 1/ζ^2
  ddz = ComplexF64(0)
  dddz = ComplexF64(0)

  for n in 1:length(C)-2
    dz -= n*C[n+2]*ζⁿ
    ζⁿ /= ζ
    dtmp = n*(n+1)*C[n+2]*ζⁿ
    ddz += dtmp
    dddz -= (n+2)*dtmp/ζ
  end
  return dz, ddz, dddz
end

function (dps::PowerSeriesDerivatives)(ζs::AbstractArray{T}) where T<:Number
  dz = zero(ζs)
  ddz = zero(ζs)
  dddz = zero(ζs)
  for (i,ζ) in enumerate(ζs)
    dz[i], ddz[i], dddz[i] = dps(ζ)
  end
  return dz, ddz, dddz
end

function power_series_first_derivative_only!(dzs::AbstractArray{T},ζs::AbstractArray{T},dps::PowerSeriesDerivatives) where {T<:Number}
  for i in eachindex(ζs)
    dzs[i] = power_series_first_derivative_only(ζs[i],dps)
  end
end

function power_series_first_derivative_only(ζs::AbstractArray{T},dps::PowerSeriesDerivatives) where {T<:Number}
  return map(ζ -> power_series_first_derivative_only(ζ,dps),ζs)
end

power_series_first_derivative_only(ζ::T,dps::PowerSeriesDerivatives) where T<:Number =
                  _power_series_first_derivative_only(ζ,dps.ps.ccoeff)

function _power_series_first_derivative_only(ζ::T,C::Vector{ComplexF64}) where T<:Number
  dz = C[1]
  ζⁿ = 1/ζ^2

  for n in 1:length(C)-2
    dz -= n*C[n+2]*ζⁿ
    ζⁿ /= ζ
  end
  return dz
end

struct PowerSeriesCache
  dz_cache :: Vector
  scale :: Vector
  dps :: PowerSeriesDerivatives
end



function evalinv_exterior(z::Vector{ComplexF64},ps::PowerSeries,
                                dps::PowerSeriesDerivatives)
#=
Evaluates the inverse of the exterior power series mapping, using a combination
of integration and Newton iteration.

This routine works well except it sometimes places points inside the
circle. This usually only happens for a flat plate.
=#

   ζ = zeros(ComplexF64,size(z))
   lenz = length(z)
   ζ0 = []
   maxiter = 10
   tol = 1.0e-8

   M = length(z)

   # sample the shape with n points, which will serve as vertices
   n = 50
   Θ = range(0,2π,length=n+1)
   prev = exp.(im.*Θ[1:n])

   # Find z values close to vertices and set ζ to the corresponding
   # prevertices
   done = zeros(Bool,size(z))

   # Now, for remaining z values, first try to integrate
   #  dζ/dt = (z - z(ζ₀))/z'(ζ) from t = 0 to t = 1,
   # with the initial condition ζ(0) = ζ₀.
   if isempty(ζ0)

     ζbase = prev

     # check for starting points on edges of the body, and rotate them
     # a bit if so
     dzbase = power_series_first_derivative_only(ζbase,dps)
     onedge = isapprox.(abs.(dzbase),0.0;atol=3*eps())
     ζbase[onedge] .*= exp(-im*1.0*(Θ[2]-Θ[1]))
     zbase = ps(ζbase)

     # Find the point on shape that is closest and use this as initial condition
     idx = []
     dist,idxtemp = findmin(abs.( repeat(transpose(z), n) - repeat(zbase, 1, M)), dims = 1)
     for k = 1:M
       push!(idx,idxtemp[k][1])
     end
     ζ0 = ζbase[idx]

     z0 = ps(ζ0)
   else
     z0 = ps(ζ0)
     if length(ζ0)==1 && lenz > 1
       ζ0 = repeat(transpose(ζ0), lenz)
       z0 = repeat(transpose(z0), lenz)
     end
     z0 = z0[.~done]
     ζ0 = ζ0[.~done]
   end
   odetol = max(tol,1e-4)
   scale = z[.~done] - z0

   dz_storage = zero(ζ0)

   p = PowerSeriesCache(dz_storage,scale,dps)

   # Checks whether an initial point will tend to make
   # the sequence of zeta values go into the circle.
   # If so, swap the initial point with its conjugate
   dζ = invfunc(ζ0,p,0.0)
   wrong_side = real(conj(ζ0).*dζ) .< 0
   ζ0[wrong_side] = conj(ζ0[wrong_side])
   z0[wrong_side] = ps(ζ0[wrong_side])
   scale = z[.~done] - z0

   p = PowerSeriesCache(dz_storage,scale,dps)

   tspan = (0.0,1.0)
   prob = ODEProblem(invfunc,ζ0,tspan,p)
   sol = solve(prob,BS3(),reltol=odetol,abstol=odetol)

   ζ[.~done] = sol.u[end] #sol.u[end][1:lenz]+im*sol.u[end][lenz+1:lenu];

   # Now use Newton iterations to improve the solution
   ζn = ζ
   k = 0
   while ~all(done) && k < maxiter
     F = z[.~done] - ps(ζn[.~done])

     ζ_notdone = ζn[.~done]
     dF = zero(ζ_notdone)
     power_series_first_derivative_only!(dF,ζ_notdone,dps)
     ζn[.~done] .= ζ_notdone .+ F./dF

     done[.~done] = abs.(F).< tol
     k += 1
   end
   F = z[.~done] - ps(ζn[.~done])
   if any(abs.(F).> tol)
     error("Check solution. Maximum residual = "*string(maximum(abs.(F))))
   end
   ζ = ζn

   return ζ

end

function invfunc(ζ,p::PowerSeriesCache,t)
  @unpack dz_cache, scale, dps = p
  power_series_first_derivative_only!(dz_cache,ζ,dps)
  return scale./dz_cache
end
