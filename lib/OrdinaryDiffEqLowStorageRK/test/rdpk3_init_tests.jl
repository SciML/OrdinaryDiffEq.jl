using Test

# Regression test for https://github.com/SciML/OrdinaryDiffEq.jl/issues/3384
#
# PIDController constructed `err = MVector{3, QT}(true, true, true)`. Because
# OrdinaryDiffEqCore imports `MVector` from `StaticArraysCore` (not
# `StaticArrays`), the Bool -> Float64 constructor is not defined and `init`
# crashes for every RDPK3 solver. The bug is masked whenever `StaticArrays` is
# loaded into the session, because `StaticArrays` adds the missing constructor
# method to `StaticArraysCore.MVector`. The existing RDPK3 tests in this file
# load `DiffEqDevTools`/`ODEProblemLibrary`, which transitively load
# `StaticArrays`, so they never exercised the failing path. The bug only
# surfaced once allocation tests were added in a leaner environment.
#
# To guard against a regression we re-run the minimal reproducer from the
# issue in a fresh Julia process that only loads `OrdinaryDiffEqLowStorageRK`.

const RDPK3_ALGS = [
    "RDPK3Sp35", "RDPK3SpFSAL35",
    "RDPK3Sp49", "RDPK3SpFSAL49",
    "RDPK3Sp510", "RDPK3SpFSAL510",
]

@testset "RDPK3 init without StaticArrays loaded (issue #3384)" begin
    for alg in RDPK3_ALGS
        script = """
        using OrdinaryDiffEqLowStorageRK
        @assert !haskey(Base.loaded_modules, Base.PkgId(Base.UUID("90137ffa-7385-5640-81b9-e52037218182"), "StaticArrays"))
        f!(du, u, p, t) = (du[1] = -0.5 * u[1]; du[2] = -1.5 * u[2])
        prob = ODEProblem(f!, [1.0, 1.0], (0.0, 1.0))
        init(prob, $alg(), dt = 0.1)
        """
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script`
        @test success(pipeline(cmd; stdout = devnull, stderr = devnull))
    end
end
