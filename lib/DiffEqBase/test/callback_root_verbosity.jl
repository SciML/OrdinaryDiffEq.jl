using DiffEqBase, SciMLBase, Test
using DiffEqBase: DEVerbosity, find_root, root_verbosity
using BracketingNonlinearSolve: NonlinearVerbosity
using SciMLLogging: None, Minimal, Standard, Detailed, All, Silent

# The root solve that locates a continuous callback's event follows the integrator's nonlinear
# verbosity (its `DEVerbosity` sub-specifier), as the integrator's other nonlinear solves do,
# instead of the bracketing solver's default logging. The root found is the same.
struct OptsOnly{O}
    opts::O
end
integ(v) = OptsOnly((; verbose = v))
@testset "callback root solve verbosity" begin
    for p in (None(), Minimal(), Standard(), Detailed(), All())
        @test root_verbosity(integ(DEVerbosity(p))) === (; verbose = NonlinearVerbosity(DEVerbosity(p).nonlinear_verbosity))
    end
    nv = NonlinearVerbosity(non_enclosing_interval = Silent())
    @test root_verbosity(integ(DEVerbosity(nonlinear_verbosity = nv))) === (; verbose = nv)
    @test root_verbosity(OptsOnly((; abstol = 1.0))) == (;)
    @test root_verbosity(nothing) == (;)
    f(t, p = nothing) = t - 0.3
    for rf in (SciMLBase.LeftRootFind, SciMLBase.RightRootFind)
        ref = find_root(f, (0.0, 1.0), rf)
        for p in (None(), Minimal(), Standard())
            @test find_root(f, (0.0, 1.0), rf; root_verbosity(integ(DEVerbosity(p)))...) === ref
        end
    end
    # a bracket that encloses no root warns at the default and Standard verbosity, and is silent at None
    g(t, p = nothing) = t + 1
    @test_logs (:warn,) match_mode = :any find_root(g, (0.0, 1.0), SciMLBase.LeftRootFind; root_verbosity(integ(DEVerbosity()))...)
    @test_logs (:warn,) match_mode = :any find_root(g, (0.0, 1.0), SciMLBase.LeftRootFind; root_verbosity(integ(DEVerbosity(Standard())))...)
    @test_logs find_root(g, (0.0, 1.0), SciMLBase.LeftRootFind; root_verbosity(integ(DEVerbosity(None())))...)
end

@testset "find_root with an explicit verbose" begin
    # the root of t - 0.3 is 0.3; ModAB with abstol = reltol = 0 narrows the bracket to adjacent floats
    f(t, p = nothing) = t - 0.3
    for rf in (SciMLBase.LeftRootFind, SciMLBase.RightRootFind), v in (None(), Standard())
        @test abs(find_root(f, (0.0, 1.0), rf; verbose = NonlinearVerbosity(v)) - 0.3) <= 2eps(0.3)
    end
    # a bracket that encloses no root warns at Standard and is silent at None
    g(t, p = nothing) = t + 1
    @test_logs (:warn,) match_mode = :any find_root(g, (0.0, 1.0), SciMLBase.LeftRootFind; verbose = NonlinearVerbosity(Standard()))
    @test_logs find_root(g, (0.0, 1.0), SciMLBase.LeftRootFind; verbose = NonlinearVerbosity(None()))
    # without a `verbose` the warning stays (the bracketing solver's default)
    @test_logs (:warn,) match_mode = :any find_root(g, (0.0, 1.0), SciMLBase.LeftRootFind)
end
