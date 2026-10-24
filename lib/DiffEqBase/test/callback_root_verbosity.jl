using DiffEqBase, SciMLBase, Test
using DiffEqBase: DEVerbosity, find_root, root_verbosity
using SciMLLogging: None, Standard

# The root solve that locates a continuous callback's event follows the integrator's nonlinear
# verbosity (its `DEVerbosity` sub-specifier), as the integrator's other nonlinear solves do,
# instead of the bracketing solver's default logging. The root found is the same.
struct OptsOnly{O}
    opts::O
end
@testset "callback root solve verbosity" begin
    @test root_verbosity(OptsOnly((; verbose = DEVerbosity(None())))) == (; verbose = DEVerbosity(None()).nonlinear_verbosity)
    @test root_verbosity(OptsOnly((; verbose = DEVerbosity(Standard())))) == (; verbose = DEVerbosity(Standard()).nonlinear_verbosity)
    @test root_verbosity(OptsOnly((; abstol = 1.0))) == (;)
    @test root_verbosity(nothing) == (;)
    f(t, p = nothing) = t - 0.3
    for rf in (SciMLBase.LeftRootFind, SciMLBase.RightRootFind)
        ref = find_root(f, (0.0, 1.0), rf)
        @test find_root(f, (0.0, 1.0), rf; root_verbosity(OptsOnly((; verbose = DEVerbosity(None()))))...) === ref
        @test find_root(f, (0.0, 1.0), rf; root_verbosity(OptsOnly((; verbose = DEVerbosity(Standard()))))...) === ref
    end
end
