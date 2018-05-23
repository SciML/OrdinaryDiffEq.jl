using DiffEqBase
using OrdinaryDiffEq
using Base.Test

@testset "Export tests" begin
  @test DiffEqBase.undefined_exports(OrdinaryDiffEq) == []
end
