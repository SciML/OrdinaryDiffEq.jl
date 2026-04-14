using Test

# Regression test for https://github.com/SciML/OrdinaryDiffEq.jl/issues/3384
#
# `PIDController` used an `MVector{3, QT}` field for its error history. Because
# `OrdinaryDiffEqCore` imports `MVector` from `StaticArraysCore` (not the full
# `StaticArrays`), neither the variadic constructor nor `setindex!` is defined
# for `MVector` in a minimal environment — both live in `StaticArrays`. So
# every RDPK3 solver crashed on `init` (construction) and would have crashed
# again on the first `step!` (setindex!) if it had gotten that far. The bug
# was masked whenever `StaticArrays` was loaded in the session, because its
# package init injects the missing methods back onto `StaticArraysCore.MVector`.
#
# The existing RDPK3 tests in `ode_low_storage_rk_tests.jl` load
# `DiffEqDevTools`/`ODEProblemLibrary`, which transitively load `StaticArrays`,
# so they never exercised the failing path. The bug only surfaced once
# allocation tests were added in a leaner environment.
#
# To guard against regressions, re-run the minimal reproducer from the issue
# in a fresh Julia process that only loads `OrdinaryDiffEqLowStorageRK`, and
# actually `solve` so stepping (not just construction) is exercised.

const RDPK3_ALGS = [
    "RDPK3Sp35", "RDPK3SpFSAL35",
    "RDPK3Sp49", "RDPK3SpFSAL49",
    "RDPK3Sp510", "RDPK3SpFSAL510",
]

@testset "RDPK3 solve without StaticArrays loaded (issue #3384)" begin
    for alg in RDPK3_ALGS
        script = """
        using OrdinaryDiffEqLowStorageRK
        @assert !haskey(Base.loaded_modules, Base.PkgId(Base.UUID("90137ffa-7385-5640-81b9-e52037218182"), "StaticArrays"))
        f!(du, u, p, t) = (du[1] = -0.5 * u[1]; du[2] = -1.5 * u[2])
        prob = ODEProblem(f!, [1.0, 1.0], (0.0, 1.0))
        sol = solve(prob, $alg(), reltol = 1.0e-6, abstol = 1.0e-8)
        sol.retcode == ReturnCode.Success || error("solve failed: \$(sol.retcode)")
        """
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script`
        @test success(pipeline(cmd; stdout = devnull, stderr = devnull))
    end
end
