using DiffEqBase, Test
using SciMLBase

# Verify that the DEVerbosity paths of `_process_verbose_param` are type-stable
# and that passing a Bool throws (v7 breaking change).

@testset "_process_verbose_param inference" begin
    # DEVerbosity passthrough is type-stable.
    @inferred DiffEqBase._process_verbose_param(DiffEqBase.DEFAULT_VERBOSE)

    # Preset dispatch is type-stable.
    @inferred DiffEqBase._process_verbose_param(DiffEqBase.SciMLLogging.Standard())

    # Bool is no longer accepted — throws ArgumentError.
    @test_throws ArgumentError DiffEqBase._process_verbose_param(true)
    @test_throws ArgumentError DiffEqBase._process_verbose_param(false)
end

# Regression test for the BoundaryValueDiffEqCore precompile failure:
# `DiffEqBase.solve`/`init` must not forward a `DEVerbosity` into kwargs
# when the caller omits the `verbose` kwarg. Downstream solver packages
# (e.g. `BoundaryValueDiffEqCore`) have their own verbosity types (e.g.
# `BVPVerbosity`) and their own default; forwarding `DEVerbosity` would
# force their `_process_verbose_param` dispatch tables to match against
# a foreign type and fail with `MethodError`.
@testset "solve/init does not force DEVerbosity into kwargs" begin
    # Fake problem + alg + __init/__solve that capture their received kwargs
    # so we can assert what DiffEqBase.solve/init actually forwards.
    Base.@kwdef mutable struct _CapturedKwargs
        solve::Union{Nothing, NamedTuple} = nothing
        init::Union{Nothing, NamedTuple} = nothing
    end

    struct _CaptureAlg <: SciMLBase.AbstractDEAlgorithm end

    struct _CaptureProb <: SciMLBase.AbstractDEProblem
        u0::Vector{Float64}
        p::SciMLBase.NullParameters
        tspan::Tuple{Float64, Float64}
        f::Nothing
        kwargs::NamedTuple
    end
    _CaptureProb() = _CaptureProb([1.0], SciMLBase.NullParameters(), (0.0, 1.0), nothing, NamedTuple())

    captured = _CapturedKwargs()

    SciMLBase.__solve(::_CaptureProb, ::_CaptureAlg; kwargs...) = (captured.solve = NamedTuple(kwargs); nothing)
    SciMLBase.__init(::_CaptureProb, ::_CaptureAlg; kwargs...) = (captured.init = NamedTuple(kwargs); nothing)
    DiffEqBase.get_concrete_problem(p::_CaptureProb, isadapt; kwargs...) = p
    DiffEqBase.isadaptive(::_CaptureAlg) = false
    DiffEqBase.prepare_alg(alg::_CaptureAlg, u0, p, prob) = alg
    DiffEqBase.check_prob_alg_pairing(::_CaptureProb, ::_CaptureAlg) = true
    SciMLBase.has_kwargs(::_CaptureProb) = false

    prob = _CaptureProb()
    alg = _CaptureAlg()

    # No explicit verbose -> must NOT be in kwargs.
    solve(prob, alg; wrap = Val(false))
    @test !haskey(captured.solve, :verbose)

    init(prob, alg)
    @test !haskey(captured.init, :verbose)

    # Explicit preset -> MUST be forwarded.
    solve(prob, alg; wrap = Val(false), verbose = DiffEqBase.SciMLLogging.Standard())
    @test haskey(captured.solve, :verbose)
    @test captured.solve.verbose isa DiffEqBase.SciMLLogging.AbstractVerbosityPreset

    # Explicit DEVerbosity -> MUST be forwarded.
    init(prob, alg; verbose = DiffEqBase.DEFAULT_VERBOSE)
    @test haskey(captured.init, :verbose)
    @test captured.init.verbose isa DiffEqBase.DEVerbosity
end
