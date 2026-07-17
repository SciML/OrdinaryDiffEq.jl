using SafeTestsets
using Pkg
using SciMLTesting

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_sparse_env()
    Pkg.activate(joinpath(@__DIR__, "sparse"))
    Pkg.develop(PackageSpec(path = dirname(dirname(dirname(@__DIR__)))))
    return Pkg.instantiate()
end

function activate_modelingtoolkit_env()
    Pkg.activate(joinpath(@__DIR__, "modelingtoolkit"))
    Pkg.develop(PackageSpec(path = dirname(dirname(dirname(@__DIR__)))))
    return Pkg.instantiate()
end

function activate_qa_env()
    return activate_group_env(joinpath(@__DIR__, "qa"); parent = [dirname(@__DIR__), joinpath(@__DIR__, "..", "..", "..")])
end

# Run QA tests (JET, Aqua)
if TEST_GROUP ∉ ("Core", "Sparse", "ModelingToolkit") && isempty(VERSION.prerelease)
    activate_qa_env()
    @time @safetestset "JET Tests" include("qa/jet.jl")
    @time @safetestset "Aqua" include("qa/qa.jl")
end

# Run functional tests
if TEST_GROUP ∉ ("QA", "Sparse", "ModelingToolkit")
    @time @safetestset "DAE jacobian2W sparse" include("dae_jacobian2w_sparse_tests.jl")
    @time @safetestset "nzval helpers" include("nzval_helpers_tests.jl")
    @time @safetestset "OOP J_t Tracking" include("oop_jt_tracking_test.jl")
    @time @safetestset "Differentiation Trait Tests" include("differentiation_traits_tests.jl")
    @time @safetestset "Autodiff Error Tests" include("autodiff_error_tests.jl")
    @time @safetestset "No Jac Tests" include("nojac_tests.jl")
    @time @safetestset "Stale W Linear Operator Tests" include("stale_w_linear_operator_tests.jl")
end

# Run sparse tests (separate environment due to ComponentArrays dep conflicts)
if TEST_GROUP == "Sparse" || TEST_GROUP == "ALL"
    activate_sparse_env()
    @time @safetestset "Non-Full Diagonal Sparsity Tests" include("sparse/nonfulldiagonal_sparse_tests.jl")
end

# Run ModelingToolkit tests (separate environment due to heavy MTK dependency)
if TEST_GROUP == "ModelingToolkit" && isempty(VERSION.prerelease)
    activate_modelingtoolkit_env()
    @time @safetestset "Jacobian Tests" include("modelingtoolkit/jacobian_tests.jl")
end
