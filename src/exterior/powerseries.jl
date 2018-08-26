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

(ps::PowerSeries)(ζ::Vector{T}) where T<:Number = ps.(ζ)

struct PowerSeriesDerivatives
  ps :: PowerSeries
end

function (dps::PowerSeriesDerivatives)(ζ::T) where T<:Number
  C = dps.ps.ccoeff
  dz = C[1]
  ζⁿ = 1/ζ^2
  ddz = ComplexF64(0)
  for n in 1:length(C)-2
    dz -= n*C[n+2]*ζⁿ
    ζⁿ /= ζ
    ddz += n*(n+1)*C[n+2]*ζⁿ
  end
  return dz, ddz
end

function (dps::PowerSeriesDerivatives)(ζs::Vector{T}) where T<:Number
  dz = zero(ζs)
  ddz = zero(ζs)
  for (i,ζ) in enumerate(ζs)
    dz[i], ddz[i] = dps(ζ)
  end
  return dz, ddz
end

function evalinv_exterior(z::Vector{ComplexF64},ps::PowerSeries,
                                dps::PowerSeriesDerivatives)
#=
Evaluates the inverse of the exterior power series mapping, using a combination
of integration and Newton iteration.
=#

   ζ = zeros(ComplexF64,size(z))
   lenz = length(z)
   ζ0 = []
   maxiter = 10
   tol = 1.0e-8

   # Find z values close to vertices and set ζ to the corresponding
   # prevertices
   done = zeros(Bool,size(z))

   # Now, for remaining z values, first try to integrate
   #  dζ/dt = (z - z(ζ₀))/z'(ζ) from t = 0 to t = 1,
   # with the initial condition ζ(0) = ζ₀.
   if isempty(ζ0)
     # choose a point on the unit circle
     ζ0 = exp.(im*zeros(lenz))
     ζ0[isapprox.(angle.(z),π;atol=eps())] .= exp(im*π)
     dz0,ddz0 = dps(ζ0)
     # check for starting points on edges of the body, and rotate them
     # a bit if so
     onedge = isapprox.(abs.(dps(ζ0)[1]),0.0;atol=eps())
     ζ0[onedge] .*= exp(im*π/20)
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
   odetol = max(tol,1e-3)
   scale = z[.~done] - z0

   ζ0 = [real(ζ0);imag(ζ0)]

   f(ζ,p,t) = invfunc(ζ,scale,dps)
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
     F = z[.~done] - ps(ζn[.~done])
     M = length(F)
     dF, ddz = dps(ζn[.~done])
     ζn[.~done] = ζn[.~done] + F./dF

     done[.~done] = abs.(F).< tol
     k += 1
   end
   F = z[.~done] - ps(ζn[.~done])
   if any(abs.(F).> tol)
     error("Check solution")
   end
   ζ = ζn

end

function invfunc(u,scale,dps::PowerSeriesDerivatives)
    lenu = length(u)
    lenzp = Int(lenu/2)
    ζ = u[1:lenzp]+im*u[lenzp+1:lenu]

    dz, ddz = dps(ζ)
    f = scale./dz
    zdot = [real(f);imag(f)]
end
