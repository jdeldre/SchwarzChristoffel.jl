@testset "Exterior map" begin
  z1 = -1.0+0.0im
  z2 = exp(im*5π/4)
  sing1 = 2
  beta = [0.5,0.5,0.5,0.5]
  z = [1.0im,-1.0,-1.0im,1.0]
  nqpts = 8
  qdat = SchwarzChristoffel.Integration.qdata(beta,nqpts)
  dquad = SchwarzChristoffel.Exterior.DQuad(z,beta,qdat)
  @test dquad(z1,z2,sing1) ≈ -0.599070117 - 0.599070117im
  dabsquad = SchwarzChristoffel.Exterior.DabsQuad(z,beta,qdat)
  @test dabsquad(z1,z2,sing1) ≈ 0.847213085

  p = SchwarzChristoffel.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
  map = SchwarzChristoffel.ExteriorMap(p)
  prev, C = SchwarzChristoffel.parameters(map)
  @test prev[2] ≈ -0.902383088-0.430934755im

  zeta = [0.1,0.5-0.75im,-0.25-0.3im]
  @test SchwarzChristoffel.evaluate(zeta,map) ≈
        Complex128[-6.934402676 - 7.689645145im,
                    0.043977374 - 1.112493651im,
                    2.411813129 - 0.044778980im]

end
