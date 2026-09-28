import OrdinaryDiffEqSDIRK
using JET

@testset "JET Tests" begin
    # JET 0.12 reports the stage variables of the unrolled out-of-place IMEX step, which are
    # assigned under `s >= n` guards, as possibly undefined.
    test_package(
        OrdinaryDiffEqSDIRK, target_modules = (OrdinaryDiffEqSDIRK,), mode = :typo,
        broken = VERSION >= v"1.12"
    )
end
