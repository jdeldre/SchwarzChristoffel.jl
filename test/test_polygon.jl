@testset "Polygon" begin
  xvert = [-1.0,1.0,0.5,-1.0]
  yvert = [-1.0,-1.0,1.0,1.0]
  p = SchwarzChristoffel.Polygon(xvert,yvert)
  @test length(SchwarzChristoffel.vertex(p)) == 4
end
