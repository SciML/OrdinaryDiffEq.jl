using OrdinaryDiffEqBDF
using ADTypes: AutoFiniteDiff
using SparseArrays
import SciMLBase
using Test

# In-place `DAEProblem` with a sparse `jac_prototype` used to `FieldError` in
# `prepare_user_sparsity` / `build_jac_config` / `sparsity_colorvec` because those
# shared ODE helpers read `f.mass_matrix` and a `DAEFunction` has no such field
# (SciML/OrdinaryDiffEq.jl#1966).

function f!(res, du, u, p, t)
    res[1] = du[1] + u[1] - u[2]
    res[2] = u[1] + u[2] - 1
    return nothing
end

const u0 = [0.0, 1.0]
const du0 = [1.0, 0.0]
const tspan = (0.0, 1.0)
const differential_vars = [true, false]

sparse_prob() = DAEProblem(
    DAEFunction(f!; jac_prototype = sparse([1.0 1.0; 1.0 1.0])),
    du0, u0, tspan; differential_vars
)

@testset "$name with sparse jac_prototype (#1966)" for (name, alg) in (
        ("DFBDF", DFBDF()),
        ("DFBDF, AutoFiniteDiff", DFBDF(autodiff = AutoFiniteDiff())),
        ("DImplicitEuler", DImplicitEuler()),
        ("DABDF2", DABDF2()),
    )
    sol = solve(sparse_prob(), alg; abstol = 1.0e-8, reltol = 1.0e-8)
    @test SciMLBase.successful_retcode(sol)
end
