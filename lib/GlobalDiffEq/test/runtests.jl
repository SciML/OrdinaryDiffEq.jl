using SafeTestsets
using Test
using SciMLTesting

run_tests(;
    env = "ODEDIFFEQ_TEST_GROUP",
    core = function ()
        @safetestset "Basic functionality" begin
            include("basic_functionality_tests.jl")
        end
        @safetestset "Algorithm traits forwarding" begin
            include("algorithm_traits_tests.jl")
        end
        @safetestset "Adjoint error estimation and control" begin
            include("adjoint_tests.jl")
        end
        @safetestset "GLEE solvers" begin
            include("glee_tests.jl")
        end
        return @safetestset "BigFloat support" begin
            include("bigfloat_tests.jl")
        end
    end,
    # QA group: Aqua + ExplicitImports + JET via SciMLTesting's run_qa. Runs in its
    # own sub-env (test/qa/Project.toml). Kept out of the curated "All" (env-bearing),
    # matching the prior JET group's exclusion.
    qa = (;
        env = joinpath(@__DIR__, "qa"),
        body = joinpath(@__DIR__, "qa", "qa.jl"),
    ),
    all = ["Core"],
)
