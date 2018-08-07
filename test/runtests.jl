using Test
using TestSetExtensions

using SchwarzChristoffel

@testset ExtendedTestSet "All tests" begin
  @includetests ARGS
end

# if isempty(ARGS)
#     include("../docs/make.jl")
# end
