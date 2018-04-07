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
        Complex128[67.206798744 + 76.6283835014im,
                   -1.116655339 + 0.5445759699im,
                    3.991294815 - 5.3064069537im]
  zeta = 0.5-0.75im
  @test SchwarzChristoffel.evaluate(zeta,map,true) ≈
        0.043977374 - 1.112493651im
  (dz,ddz) = SchwarzChristoffel.evalderiv(zeta,map,true)
  @test dz ≈ -1.116655339 + 0.5445759699im
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
  zeta = 1.0+3.0im
  @test SchwarzChristoffel.evaluate(zeta,map,false) ≈
               0.816139620 + 3.029559043im
  dz, ddz = SchwarzChristoffel.evalderiv(zeta,map,false)
  @test dz ≈ 1.030528002 + 0.004449924im

  ccoeff, dcoeff = SchwarzChristoffel.coefficients(map)
  @test ccoeff[1] ≈ 1.019795410

end

@testset "Power series map" begin

  a = 1
  b = 1
  c = Complex128[0.5(a+b),0,0.5(a-b)]
  m = SchwarzChristoffel.PowerMap(c)
  @test area(m) ≈    3.141592654
  @test isapprox(centroid(m),Complex128(0);atol=eps())
  @test Jmoment(m) ≈ 1.570796327
  a = 1
  b = 0
  c = Complex128[0.5(a+b),0,0.5(a-b)]
  m = SchwarzChristoffel.PowerMap(c)
  @test area(m) ≈    Float64(0)
  @test isapprox(centroid(m),Complex128(0);atol=eps())
  @test Jmoment(m) ≈ Float64(0)

  m = SchwarzChristoffel.PowerMap(Complex128[1,0,0,0,0,0.1])
  @test area(m) ≈ 3.015928947
  @test isapprox(centroid(m),Complex128(0);atol=eps())
  @test Jmoment(m) ≈ 1.475920229

end
