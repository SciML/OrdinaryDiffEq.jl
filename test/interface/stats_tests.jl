# stats.nf tests
using OrdinaryDiffEq, Test, ADTypes
x = Ref(0)
function f(u, p, t)
    x[] += 1
    return 5 * u
end
function g(du, u, p, t)
    x[] += 1
    @. du = 5 * u
end

u0 = [1.0, 1.0]
tspan = (0.0, 1.0)
probop = ODEProblem(f, u0, tspan)
probip = ODEProblem(g, u0, tspan)

@testset "stats_tests" begin
    @testset "$prob" for prob in [probop, probip]
        @testset "$alg" for alg in [BS3, Tsit5, Vern7, Vern9, ROCK4]
            x[] = 0
            sol = solve(prob, alg())
            @test x[] == sol.stats.nf
        end
        @testset "$alg" for alg in [Rodas5P, KenCarp4]
            @testset "$kwargs" for kwargs in [(autodiff = AutoForwardDiff(),),
                (autodiff = AutoFiniteDiff(fdtype = Val{:forward}()),),
                (autodiff = AutoFiniteDiff(fdtype = Val{:central}()),),
                (autodiff = AutoFiniteDiff(fdtype = Val{:complex}()),)]
                x[] = 0
                sol = solve(prob, alg(; kwargs...))
                @test x[] == sol.stats.nf
            end
        end
    end
end
