using DiffEqBase
import SciMLBase
using Test

@test DiffEqBase.undefined_exports(DiffEqBase) == []

@testset "Problem type bindings" begin
    for name in (
            :AbstractLinearProblem, :AbstractNonlinearProblem,
            :AbstractIntegralProblem, :AbstractSteadyStateProblem,
        )
        @testset "$name" begin
            expected = getproperty(SciMLBase, name)
            @test getproperty(DiffEqBase, name) === expected
            @test !Base.isexported(DiffEqBase, name)
            consumer = Module(gensym(:DiffEqBaseConsumer))
            Core.eval(consumer, :(import DiffEqBase: $name))
            @test Core.eval(consumer, name) === expected
            Core.eval(consumer, :(accepts(::$name) = true))
            @test hasmethod(Core.eval(consumer, :accepts), Tuple{expected})
        end
    end
end
