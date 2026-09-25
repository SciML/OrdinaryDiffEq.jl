using OrdinaryDiffEq
using Test

import OrdinaryDiffEqCore
import OrdinaryDiffEqExtrapolation
import OrdinaryDiffEqFIRK
import OrdinaryDiffEqPDIRK
import OrdinaryDiffEqPRK

const THREADING_OPTIONS = (:Sequential, :BaseThreads, :PolyesterThreads)
const THREADING_MODULES = (
    OrdinaryDiffEq, OrdinaryDiffEqExtrapolation, OrdinaryDiffEqFIRK,
    OrdinaryDiffEqPDIRK, OrdinaryDiffEqPRK,
)

@testset "threading options resolve qualified from every module taking `threading`" begin
    for mod in THREADING_MODULES, name in THREADING_OPTIONS
        @test isdefined(mod, name)
        # Same binding as OrdinaryDiffEqCore's, not a shadowing copy.
        @test getproperty(mod, name) === getproperty(OrdinaryDiffEqCore, name)
    end

    # v6 reached these as `OrdinaryDiffEq.PolyesterThreads()`; that is the spelling
    # v7 broke by dropping the umbrella's import, and the one being restored.
    @test OrdinaryDiffEq.PolyesterThreads() isa OrdinaryDiffEqCore.AbstractThreadingOption
    @test OrdinaryDiffEqExtrapolation.PolyesterThreads() === OrdinaryDiffEq.PolyesterThreads()
end

@static if VERSION >= v"1.11.0-DEV.469"
    @testset "threading options are declared public, not exported" begin
        for mod in THREADING_MODULES, name in THREADING_OPTIONS
            @test Base.ispublic(mod, name)
            # `public`, deliberately not `export`: `Sequential` is too generic to put
            # in every user's namespace, and bare use never worked in v6 either.
            # (`names` lists public names too, so `isexported` is the right check.)
            @test !Base.isexported(mod, name)
        end
    end
end

@testset "bare use stays an error" begin
    # This file's `using OrdinaryDiffEq` would bind these if they were exported.
    for name in THREADING_OPTIONS
        @test !isdefined(@__MODULE__, name)
    end
end

@testset "the qualified options drive the solvers" begin
    prob = ODEProblem((du, u, p, t) -> (du[1] = -u[1]; nothing), [1.0], (0.0, 1.0))
    for threading in (
            false, OrdinaryDiffEq.Sequential(), OrdinaryDiffEq.BaseThreads(),
        )
        alg = OrdinaryDiffEqExtrapolation.ExtrapolationMidpointDeuflhard(; threading)
        @test successful_retcode(solve(prob, alg))
    end
end
