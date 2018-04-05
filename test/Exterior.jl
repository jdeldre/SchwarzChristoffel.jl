@testset "Exterior map" begin
  z1 = -1.0+0.0im
  z2 = exp(im*5π/4)
  sing1 = 2
  beta = [0.5,0.5,0.5,0.5]
  z = [1.0im,-1.0,-1.0im,1.0]
  nqpts = 8
  qdat = SchwarzChristoffel.Integration.qdata(beta,nqpts)
  dquad = SchwarzChristoffel.Exterior.DQuad(beta,qdat)
  @test dquad([z1],[z2],[sing1],z)[1] ≈ -0.599070117 - 0.599070117im
  dabsquad = SchwarzChristoffel.Exterior.DabsQuad(beta,qdat)
  @test dabsquad([z1],[z2],[sing1],z)[1] ≈ 0.847213085

  p = SchwarzChristoffel.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
  map = SchwarzChristoffel.ExteriorMap(p)
  prev, C = SchwarzChristoffel.parameters(map)
  @test prev[3] ≈ -0.902383088-0.430934755im

  zeta = [0.1,0.5-0.75im,-0.25-0.3im]
  @test SchwarzChristoffel.evaluate(zeta,map,true) ≈
        Complex128[-6.934402676 - 7.689645145im,
                    0.043977374 - 1.112493651im,
                    2.411813129 - 0.044778980im]
  (dz,ddz) = SchwarzChristoffel.evalderiv(zeta,map,true)
  @test dz ≈
        Complex128[67.206798225 + 76.628383817im,
                   -1.116655325 + 0.544575954im,
                    3.991294846 - 5.306406913im]
  zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im]
  @test SchwarzChristoffel.evaluate(zeta,map,false) ≈
        Complex128[ 0.816139620 + 3.029559043im,
                   -2.252366325 - 2.085230469im,
                   -0.333104076 + 0.975837123im]
  (dz,ddz) = SchwarzChristoffel.evalderiv(zeta,map,false)
  @test dz ≈
        Complex128[ 1.030528002 + 0.004449924im,
                    1.006955879 - 0.011501136im,
                    1.300781073 - 0.266624652im]

  ccoeff, dcoeff = SchwarzChristoffel.coefficients(map)
  @test ccoeff[1] ≈ 1.019795410

end
