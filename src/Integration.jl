module Integration

using Compat
using FastGaussQuadrature

export qdata

function qdata(beta::Vector{Float64},nqpts::Int64)

  n = length(beta)
  qnode = zeros(nqpts,n+1)
  qwght = zeros(nqpts,n+1)
  for j = findall(beta.>-1)
    qnode[:,j],qwght[:,j] = gaussjacobi(nqpts,0,beta[j])
  end
  qnode[:,n+1],qwght[:,n+1] = gaussjacobi(nqpts,0,0)
  return qnode, qwght


end


end
