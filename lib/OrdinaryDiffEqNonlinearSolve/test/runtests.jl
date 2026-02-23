using SafeTestsets
using Pkg

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_modelingtoolkit_env()
    Pkg.activate(joinpath(@__DIR__, "modelingtoolkit"))
    lib_dir = dirname(dirname(@__DIR__))
    Pkg.develop(
        [
            PackageSpec(path = dirname(dirname(dirname(@__DIR__)))),
            PackageSpec(path = dirname(@__DIR__)),
            PackageSpec(path = joinpath(lib_dir, "OrdinaryDiffEqBDF")),
            PackageSpec(path = joinpath(lib_dir, "OrdinaryDiffEqSDIRK")),
        ]
    )
    return Pkg.instantiate()
end

# Run functional tests
if TEST_GROUP ∉ ("QA", "ModelingToolkit")
    @time @safetestset "Newton Tests" include("newton_tests.jl")
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
if TEST_GROUP ∉ ("Core", "ModelingToolkit") && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end

# Run ModelingToolkit tests (separate environment due to heavy MTK dependency)
if TEST_GROUP == "ModelingToolkit" && isempty(VERSION.prerelease)
    activate_modelingtoolkit_env()
    @time @safetestset "NLStep Tests" include("modelingtoolkit/nlstep_tests.jl")
    @time @safetestset "Preconditioner Tests" include("modelingtoolkit/preconditioners.jl")
    @time @safetestset "DAE Initialize Integration" include("modelingtoolkit/dae_initialize_integration.jl")
end
