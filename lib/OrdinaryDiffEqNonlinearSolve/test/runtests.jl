using SafeTestsets

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

# Run functional tests
if TEST_GROUP != "QA"
    @time @safetestset "Newton Tests" include("newton_tests.jl")
    @time @safetestset "Sparse Algebraic Detection" include("sparse_algebraic_detection_tests.jl")
    @time @safetestset "Sparse DAE Initialization" include("sparse_dae_initialization_tests.jl")
    @time @safetestset "Linear Nonlinear Solver Tests" include("linear_nonlinear_tests.jl")
    @time @safetestset "Linear Solver Tests" include("linear_solver_tests.jl")
    @time @safetestset "Linear Solver Split ODE Tests" include("linear_solver_split_ode_tests.jl")
    @time @safetestset "Mass Matrix Tests" include("mass_matrix_tests.jl")
    @time @safetestset "W-Operator Prototype Tests" include("wprototype_tests.jl")
    @time @safetestset "DAE Initialization Tests" include("dae_initialization_tests.jl")
    @time @safetestset "CheckInit Tests" include("checkinit_tests.jl")
end

# Run QA tests (JET, Aqua)
if TEST_GROUP != "Core" && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end
