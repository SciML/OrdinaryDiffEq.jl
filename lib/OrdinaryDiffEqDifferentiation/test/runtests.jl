using SafeTestsets
using Pkg

const TEST_GROUP = get(ENV, "ODEDIFFEQ_TEST_GROUP", "ALL")

function activate_sparse_env()
    Pkg.activate(joinpath(@__DIR__, "sparse"))
    # Develop the top-level OrdinaryDiffEq package (which pulls all subpackages)
    Pkg.develop(PackageSpec(path = dirname(dirname(dirname(@__DIR__)))))
    return Pkg.instantiate()
end

# Run QA tests (JET, Aqua)
if TEST_GROUP ∉ ("Core", "Sparse") && isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
end

# Run functional tests
if TEST_GROUP ∉ ("QA", "Sparse")
    @time @safetestset "OOP J_t Tracking" include("oop_jt_tracking_test.jl")
    @time @safetestset "Differentiation Trait Tests" include("differentiation_traits_tests.jl")
    @time @safetestset "Autodiff Error Tests" include("autodiff_error_tests.jl")
    @time @safetestset "No Jac Tests" include("nojac_tests.jl")
end

# Run sparse tests (separate environment due to SparseConnectivityTracer/ComponentArrays dep conflicts)
if TEST_GROUP == "Sparse" || TEST_GROUP == "ALL"
    activate_sparse_env()
    @time @safetestset "Non-Full Diagonal Sparsity Tests" include("sparse/nonfulldiagonal_sparse_tests.jl")
end
