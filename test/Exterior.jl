@testset "Exterior map" begin
  z1 = -1.0+0.0im
  z2 = exp(im*5π/4)
  sing1 = 2
  beta = [0.5,0.5,0.5,0.5]
  z = [1.0im,-1.0,-1.0im,1.0]
  nqpts = 8
  qdat = SchwarzChristoffel.Exterior.Integration.qdata(beta,nqpts)
  dquad = SchwarzChristoffel.Exterior.DQuad(beta,qdat)
  @test dquad([z1],[z2],[sing1],z)[1] ≈ -0.599070117 - 0.599070117im
  dabsquad = SchwarzChristoffel.Exterior.DabsQuad(beta,qdat)
  @test dabsquad([z1],[z2],[sing1],z)[1] ≈ 0.847213085

  p = SchwarzChristoffel.Polygon([-1.0,0.2,1.0,-1.0],[-1.0,-1.0,0.5,1.0])
  m = SchwarzChristoffel.ExteriorMap(p)
  dm = SchwarzChristoffel.DerivativeMap(m)
  m⁻¹ = SchwarzChristoffel.InverseMap(m)
  prev, C = SchwarzChristoffel.parameters(m)
  @test prev[3] ≈ -0.902383088-0.430934755im

  zeta = [0.1,0.5-0.75im,-0.25-0.3im]
  z = m(zeta;inside=true)
  @test z ≈
        Complex128[-6.934402676 - 7.689645145im,
                    0.043977374 - 1.112493651im,
                    2.411813129 - 0.044778980im]
  @test isapprox(m⁻¹(z;inside=true),zeta;atol=eps())

  (dz,ddz) = dm(zeta;inside=true)
  @test isapprox(dz,
        Complex128[67.206798744 + 76.628383501im,
                   -1.116655339 + 0.5445759699im,
                    3.991294815 - 5.3064069537im];atol=eps())
  zeta = 0.5-0.75im
  z = m(zeta;inside=true)
  @test z ≈ 0.043977374 - 1.112493651im
  @test m⁻¹(z;inside=true) ≈ zeta
  (dz,ddz) = dm(zeta;inside=true)
  @test dz ≈ -1.116655339 + 0.5445759699im
  zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im]
  z = m(zeta)
  @test z ≈
        Complex128[ 0.816139620 + 3.029559043im,
                   -2.252366325 - 2.085230469im,
                   -0.333104076 + 0.975837123im]
  @test m⁻¹(z) ≈ zeta
  (dz,ddz) = dm(zeta)
  @test dz ≈
        Complex128[ 1.030528002 + 0.004449924im,
                    1.006955879 - 0.011501136im,
                    1.300781073 - 0.266624652im]
  zeta = 1.0+3.0im
  z = m(zeta)
  @test z ≈ 0.816139620 + 3.029559043im
  @test m⁻¹(z) ≈ zeta
  dz, ddz = dm(zeta)
  @test dz ≈ 1.030528002 + 0.004449924im

  ccoeff, dcoeff = SchwarzChristoffel.coefficients(m)
  @test ccoeff[1] ≈ 1.019789030 - 4.043787137e-8im

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

  a = 1
  b = 0.1
  c = Complex128[0.5(a+b),0,0.5(a-b)]
  m = SchwarzChristoffel.PowerMap(c)
  m⁻¹ = SchwarzChristoffel.InverseMap(m)
  dm = SchwarzChristoffel.DerivativeMap(m)
  zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im]
  z = m(zeta)
  @test z ≈
      Complex128[0.595+1.515im,-1.2125-0.9875im,
                 0.195909090909im]
  @test isapprox(m⁻¹(z),zeta;atol=eps())
  dz, ddz = dm(zeta)
  @test dz ≈ Complex128[0.586 + 0.027im,
                        0.55 + 0.05625im,
                        0.921900826]
  zeta = 1.0+3.0im
  z = m(zeta)
  @test z ≈ 0.595+1.515im
  @test isapprox(m⁻¹(z),zeta;atol=eps())
  dz, ddz = dm(zeta)
  @test dz ≈ 0.586 + 0.027im

  a = 1
  b = 0
  c = Complex128[0.5(a+b),0,0.5(a-b)]
  m = SchwarzChristoffel.PowerMap(c)
  m⁻¹ = SchwarzChristoffel.InverseMap(m)
  dm = SchwarzChristoffel.DerivativeMap(m)
  z = [-1.0+0im,-1.1+0im,1.0+0im,1.1+0im]
  @test isapprox(m(m⁻¹(z)),z;atol=eps())

end

@testset "Karman-Trefftz map" begin
  ν = 1.94; ϵ = 0.08*sqrt(2); δ = 3π/4; C = 0.25;
  m = SchwarzChristoffel.KarmanTrefftzMap(ν,ϵ,δ,C)
  @test area(m) ≈ 0.08510380847113
  @test centroid(m) ≈ -0.0796138968044 + 0.0271627776676im
  @test Jmoment(m) ≈ 0.00538941353392

  zeta = [1.0+3.0im,-2.0-2.0im,0.0+1.1im]
  z = m(zeta)
  @test isapprox(z,Complex128[0.269784+0.768849im,
                      -0.616452-0.470251im,
                      -0.0310697+0.138951im];atol=1e-5)

  m⁻¹ = SchwarzChristoffel.InverseMap(m)
  @test isapprox(m⁻¹(m(zeta)),zeta;atol=eps())


end
