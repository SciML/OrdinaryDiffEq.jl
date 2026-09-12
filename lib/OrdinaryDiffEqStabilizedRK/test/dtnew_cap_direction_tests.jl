using OrdinaryDiffEqStabilizedRK, OrdinaryDiffEqCore, Test

# `dtnew_modification` caps the controller's proposed step at the method's
# stability limit. The limit is a magnitude, so `min(dtnew, bound)` silently does
# nothing for a backward (tdir < 0) solve, where dtnew < 0 < bound: the cap was
# applied going forwards and dropped going backwards.
#
# NOTE: these methods' caps are currently unreachable from `calc_dt_propose!`
# because `OrdinaryDiffEqStabilizedRK` defines `dtnew_modification` without
# importing it from `OrdinaryDiffEqCore`, so the generic identity method is what
# actually runs (see the issue linked from the PR). The tests below therefore
# call the package's own method directly.

f(u, p, t) = -1000.0 .* u

const CAPPED_ALGS = (
    "ROCK2" => ROCK2(), "ROCK4" => ROCK4(), "SERK2" => SERK2(),
    "ESERK4" => ESERK4(), "ESERK5" => ESERK5(), "RKMC2" => RKMC2(),
    "RKL1" => RKL1(), "RKL2" => RKL2(), "RKG1" => RKG1(), "RKG2" => RKG2(),
)

@testset "stability cap is direction-agnostic" begin
    for (algname, alg) in CAPPED_ALGS
        @testset "$algname" begin
            capped = map(((0.0, 1.0), (1.0, 0.0))) do tspan
                integ = init(
                    ODEProblem(f, [1.0], tspan), alg;
                    abstol = 1.0e-6, reltol = 1.0e-6
                )
                integ.eigen_est = 1.0e10 # make the stability bound bind for every method
                dtnew = integ.tdir * 1.0 # far above it
                return abs(
                    OrdinaryDiffEqStabilizedRK.dtnew_modification(integ, alg, dtnew)
                )
            end
            @test capped[1] < 1.0        # the cap actually binds forward
            @test capped[1] == capped[2] # and identically backward
        end
    end
end
