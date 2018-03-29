@testset "Polygon" begin
  xvert = [-1.0,1.0,0.5,-1.0]
  yvert = [-1.0,-1.0,1.0,1.0]
  p = SchwarzChristoffel.Polygon(xvert,yvert)
  @test length(SchwarzChristoffel.vertex(p)) == 4
  z = Complex128(0.2+im*0)
  @test SchwarzChristoffel.isinpoly(z,p)
  z = Complex128(-1+im*0)
  @test SchwarzChristoffel.isinpoly(z,p,1e-8)
  beta = SchwarzChristoffel.interiorangle(p)
  @test beta[1] ≈ 0.5
end
