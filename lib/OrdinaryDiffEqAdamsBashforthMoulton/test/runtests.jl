using SafeTestsets

@time @safetestset "ABM Convergence Tests" include("abm_convergence_tests.jl")
@time @safetestset "Adams Variable Coefficients Tests" include("adams_tests.jl")
@time @safetestset "Regression test for threading versions vs non threading versions" include("regression_test_threading.jl")
@time @safetestset "JET Tests" include("jet.jl")
@time @safetestset "Aqua" include("qa_tests.jl")
