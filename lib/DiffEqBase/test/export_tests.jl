using DiffEqBase
import SciMLBase
using Test

@test DiffEqBase.undefined_exports(DiffEqBase) == []

@testset "Qualified SciMLBase problem types" begin
    for name in (
            :AbstractLinearProblem, :AbstractNonlinearProblem,
            :AbstractSteadyStateProblem, :AbstractIntegralProblem,
        )
        @test isdefined(DiffEqBase, name)
        @test isdefined(DiffEqBase, name) && getfield(DiffEqBase, name) === getfield(SciMLBase, name)
    end
end
