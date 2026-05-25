module EnzymeAccumTangentTests

using Test
using DiffEqBase
import ChainRulesCore, Enzyme  # triggers DiffEqBaseEnzymeExt
import SciMLStructures
import SciMLStructures: Tunable

mutable struct MockMTKParams
    tunable::Vector{Float64}
    caches::Tuple{Vector{Float64}}
end

SciMLStructures.isscimlstructure(::MockMTKParams) = true
SciMLStructures.ismutablescimlstructure(::MockMTKParams) = true
function SciMLStructures.canonicalize(::Tunable, p::MockMTKParams)
    return p.tunable, (val) -> MockMTKParams(collect(val), p.caches), true
end
function SciMLStructures.replace!(::Tunable, p::MockMTKParams, val)
    p.tunable .= val
    return nothing
end

const EXT = Base.get_extension(DiffEqBase, :DiffEqBaseEnzymeExt)

@testset "EnzymeExt._accum_tangent! gates non-Tunable walk on diff_tunables (NS#936)" begin
    # Regression for SciML/NonlinearSolve.jl#936 (ported to DiffEqBase). The
    # reverse rule mirrors `sensealg.diff_tunables` into this kwarg:
    #
    #   * `diff_tunables = true`  (default) — backpass returned a
    #     Tunable-only cotangent; leave non-Tunable fields of any
    #     structured `darg` alone.
    #   * `diff_tunables = false` — backpass returned a *structured*
    #     cotangent (e.g. caches from SCC `explicitfuns!` coupling, or any
    #     MTKParameters-bound ODE with non-empty caches). Walk every
    #     non-Tunable field of `darg` into `dval` so the meaningful
    #     contribution lands.

    # diff_tunables = false: caches must accumulate.
    dval = MockMTKParams([0.0, 0.0], (zeros(3),))
    darg = MockMTKParams([1.0, 2.0], ([10.0, 20.0, 30.0],))
    EXT._accum_tangent!(dval, darg; diff_tunables = false)
    @test dval.tunable == [1.0, 2.0]
    @test dval.caches[1] == [10.0, 20.0, 30.0]

    # And `+=` semantics on repeated calls.
    darg2 = MockMTKParams([0.5, 0.5], ([1.0, 2.0, 3.0],))
    EXT._accum_tangent!(dval, darg2; diff_tunables = false)
    @test dval.tunable == [1.5, 2.5]
    @test dval.caches[1] == [11.0, 22.0, 33.0]

    # diff_tunables = true (default): only Tunable touched.
    dval2 = MockMTKParams([0.0, 0.0], (zeros(3),))
    darg3 = MockMTKParams([1.0, 2.0], ([10.0, 20.0, 30.0],))
    EXT._accum_tangent!(dval2, darg3)
    @test dval2.tunable == [1.0, 2.0]
    @test dval2.caches[1] == zeros(3)
end

end
