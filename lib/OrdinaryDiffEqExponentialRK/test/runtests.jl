using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "Linear-Nonlinear Krylov Methods Tests" include("linear_nonlinear_krylov_tests.jl")
    @time @safetestset "Linear-Nonlinear Convergence Tests" include("linear_nonlinear_convergence_tests.jl")
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "FUNCTIONAL" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
