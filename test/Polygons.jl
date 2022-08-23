@testset "Polygons" begin
  xvert = [-1.0,1.0,0.5,-1.0]
  yvert = [-1.0,-1.0,1.0,1.0]
  p = SchwarzChristoffel.Polygon(xvert,yvert)
  @test length(SchwarzChristoffel.vertex(p)) == 4
  z = ComplexF64(0.2+im*0)
  @test SchwarzChristoffel.isinpoly(z,p)
  z = ComplexF64(-1+im*0)
  @test SchwarzChristoffel.isinpoly(z,p,1e-8)
  beta = SchwarzChristoffel.interiorangle(p)
  @test beta[1] ≈ 0.5
end

@testset "Linked lines" begin

  x = -0.5:0.1:0.5
  p = SchwarzChristoffel.LinkedLines(x,0.2*sin.(4x))

  @test p.angle[1] ≈ p.angle[11] ≈ 0.0
  @test p.angle[6] ≈ p.angle[16] ≈ 1.0


end

@testset "NACA airfoil" begin

  w = SchwarzChristoffel.naca4(0.04,0.4,0.12;len=1)
  p = SchwarzChristoffel.Polygon(w)
  @test p.angle[1] ≈ 0.1130140167872078

end
