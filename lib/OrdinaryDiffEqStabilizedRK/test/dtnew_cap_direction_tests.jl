using OrdinaryDiffEqStabilizedRK, OrdinaryDiffEqCore, Test

# `dtnew_modification` caps the controller's proposed step at the method's
# stability limit. The limit is a magnitude, so `min(dtnew, bound)` silently does
# nothing for a backward (tdir < 0) solve, where dtnew < 0 < bound: the cap was
# applied going forwards and dropped going backwards.

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

# `calc_dt_propose!` dispatches on `OrdinaryDiffEqCore.dtnew_modification`. These
# methods are only reached if this package *extends* that function rather than
# defining a same-named one of its own, which is easy to get wrong by leaving it
# out of the `import` list - and silently disables every cap above.
@testset "stability cap is reachable from OrdinaryDiffEqCore" begin
    @test OrdinaryDiffEqCore.dtnew_modification ===
        OrdinaryDiffEqStabilizedRK.dtnew_modification
    for (algname, alg) in CAPPED_ALGS
        @testset "$algname" begin
            @test OrdinaryDiffEqCore.has_dtnew_modification(alg)
            integ = init(
                ODEProblem(f, [1.0], (0.0, 1.0)), alg;
                abstol = 1.0e-6, reltol = 1.0e-6
            )
            integ.eigen_est = 1.0e10
            @test abs(OrdinaryDiffEqCore.dtnew_modification(integ, alg, 1.0)) < 1.0
        end
    end
end
