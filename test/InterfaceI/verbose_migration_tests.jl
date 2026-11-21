using OrdinaryDiffEq
using Test

prob = ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), [1.0], (0.0, 1.0))
alg = Tsit5()

@testset "verbose migration path is reachable from `using OrdinaryDiffEq`" begin
    @test :DEVerbosity in names(OrdinaryDiffEq)
    @test :SciMLLogging in names(OrdinaryDiffEq)

    @test successful_retcode(solve(prob, alg; verbose = DEVerbosity(SciMLLogging.None())))
    @test successful_retcode(solve(prob, alg; verbose = DEVerbosity()))

    err = try
        solve(prob, alg; verbose = false)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError

    # `ODEVerbosity` never existed; the message must not send users after it.
    @test !occursin("ODEVerbosity", err.msg)

    # Every `solve(...)` line the message suggests must run verbatim in a scope whose
    # only import is `using OrdinaryDiffEq`, which is what this file provides.
    suggestions = filter(
        line -> occursin("verbose = DEVerbosity", line), split(err.msg, '\n')
    )
    @test length(suggestions) == 2
    for line in suggestions
        @test successful_retcode(eval(Meta.parse(strip(line))))
    end
end
