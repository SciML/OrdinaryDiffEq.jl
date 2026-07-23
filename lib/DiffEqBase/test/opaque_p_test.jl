using DiffEqBase
using SciMLBase
using RespecializeParams
using Test
using ForwardDiff
using SymbolicIndexingInterface: SymbolCache

# An `isbits` parameter struct. `get_concrete_problem` on an AutoDePSpecialize
# in-place ODEFunction with such a `p` should route through the opaque path:
# the returned problem has `prob.p :: OpaqueParams` and `prob.f.f` is a
# FunctionWrappersWrapper whose signature carries `OpaqueParams` in the p
# slot.

struct LinearP
    k::Float64
end

struct VecP
    ks::Vector{Float64}
end

function linear_rhs!(du, u, p::LinearP, t)
    @inbounds du[1] = -p.k * u[1]
    return nothing
end

deprob(u0, p) = ODEProblem{true, DiffEqBase.AutoDePSpecialize}(
    linear_rhs!, u0, (0.0, 1.0), p,
)

# Trigger the same code path the high-level `solve(...)` does, without
# requiring a concrete integrator package as a test dep.
concretize(prob, alg = nothing) = DiffEqBase.get_concrete_problem(prob, true; alg = alg)

@testset "opaque-p hook" begin
    @testset "AutoDePSpecialize opaque-ifies isbits p" begin
        cp = concretize(deprob([1.0], LinearP(0.5)))
        @test cp.p isa RespecializeParams.OpaqueParams
        # The unpack-then-call path works on the wrapped rhs:
        du = [0.0]
        cp.f(du, [1.0], cp.p, 0.0)
        @test du[1] ≈ -0.5
    end

    @testset "AutoSpecialize leaves prob.p untouched" begin
        prob = ODEProblem{true, SciMLBase.AutoSpecialize}(
            linear_rhs!, [1.0], (0.0, 1.0), LinearP(0.5),
        )
        cp = concretize(prob)
        @test cp.p isa LinearP
        # Still function-wrapped, with the concrete p type in the signature:
        @test cp.f.f isa DiffEqBase.FunctionWrappersWrappers.FunctionWrappersWrapper
    end

    @testset "FullSpecialize leaves prob.p untouched" begin
        prob = ODEProblem{true, SciMLBase.FullSpecialize}(
            linear_rhs!, [1.0], (0.0, 1.0), LinearP(0.5),
        )
        cp = concretize(prob)
        @test cp.p isa LinearP
    end

    @testset "NullParameters skipped" begin
        f_null!(du, u, p, t) = (du[1] = -u[1]; nothing)
        prob = ODEProblem{true, DiffEqBase.AutoDePSpecialize}(
            f_null!, [1.0], (0.0, 1.0),
        )
        cp = concretize(prob)
        @test cp.p isa SciMLBase.NullParameters
    end

    @testset "non-isbits p is de-specialized via OpaqueRef" begin
        # A non-isbits struct (Vector field). It packs by reference into an
        # OpaqueRef, and the wrapper signature carries OpaqueRef in the p slot.
        f_vecp!(du, u, p::VecP, t) = (@inbounds du[1] = -p.ks[1] * u[1]; nothing)
        prob = ODEProblem{true, DiffEqBase.AutoDePSpecialize}(f_vecp!, [1.0], (0.0, 1.0), VecP([0.5]))
        cp = concretize(prob)
        @test cp.p isa RespecializeParams.OpaqueRef
        du = [0.0]
        cp.f(du, [1.0], cp.p, 0.0)
        @test du[1] ≈ -0.5
        @test RespecializeParams.unpack(cp.p, VecP).ks == [0.5]

        # a bare Vector parameter is likewise de-specialized into an OpaqueRef
        f_vec!(du, u, p::Vector{Float64}, t) = (@inbounds du[1] = -p[1] * u[1]; nothing)
        cpv = concretize(ODEProblem{true, DiffEqBase.AutoDePSpecialize}(f_vec!, [1.0], (0.0, 1.0), [0.5]))
        @test cpv.p isa RespecializeParams.OpaqueRef
    end

    @testset "symbolic-system problems are declined (MTK safety)" begin
        # A problem whose `f` carries a symbolic system (`has_sys`, as every
        # ModelingToolkit problem does) must NOT be opaque-ified: its `p` has to
        # stay concrete for the initialization pipeline and symbolic parameter
        # indexing. `SymbolCache` stands in for an MTK `System` here without the
        # ModelingToolkit dependency. The opaque path declines and falls back.
        f_sym!(du, u, p, t) = (@inbounds du[1] = -u[1]; nothing)
        ff = ODEFunction(f_sym!; sys = SymbolCache([:x], [:k], :t))
        @test SciMLBase.has_sys(ff)

        # non-isbits p that would otherwise pack into an OpaqueRef:
        cp = concretize(ODEProblem{true, DiffEqBase.AutoDePSpecialize}(ff, [1.0], (0.0, 1.0), VecP([0.5])))
        @test cp.p isa VecP
        @test !(cp.p isa RespecializeParams.OpaqueRef)

        # isbits p that would otherwise pack into an OpaqueParams:
        cpi = concretize(ODEProblem{true, DiffEqBase.AutoDePSpecialize}(ff, [1.0], (0.0, 1.0), LinearP(0.5)))
        @test cpi.p isa LinearP
        @test !(cpi.p isa RespecializeParams.OpaqueParams)
    end

    @testset "already-packed p is not re-wrapped (idempotent)" begin
        cp = concretize(deprob([1.0], LinearP(0.5)))
        cp2 = concretize(cp)   # re-concretize the opaque problem
        @test cp2.p isa RespecializeParams.OpaqueParams
        @test typeof(cp2.f) === typeof(cp.f)   # no nested OpaqueVoid wrapping
    end

    @testset "two LinearP problems share the wrapped-f type" begin
        cp_a = concretize(deprob([1.0], LinearP(0.5)))
        cp_b = concretize(deprob([2.0], LinearP(1.5)))
        @test typeof(cp_a.f) === typeof(cp_b.f)
        @test typeof(cp_a.p) === typeof(cp_b.p) === RespecializeParams.OpaqueParams
    end

    @testset "problems with different isbits p types share the wrapped-f type" begin
        f_nt!(du, u, p, t) = (@inbounds du[1] = -p.k * u[1]; nothing)
        prob_nt = ODEProblem{true, DiffEqBase.AutoDePSpecialize}(
            f_nt!, [1.0], (0.0, 1.0), (k = 0.5,),
        )
        cp_struct = concretize(deprob([1.0], LinearP(0.5)))
        cp_nt = concretize(prob_nt)
        @test typeof(cp_struct.f) === typeof(cp_nt.f)
        @test typeof(cp_struct.p) === typeof(cp_nt.p) === RespecializeParams.OpaqueParams
    end

    @testset "unsafe_unpack of the packed p yields the original value" begin
        p = LinearP(0.5)
        cp = concretize(deprob([1.0], p))
        @test RespecializeParams.unsafe_unpack(cp.p, LinearP) === p
    end

    # The whole point of AutoDePSpecialize is that the de-specialized RHS still
    # runs a type-stable, allocation-free inner call: the FunctionWrappersWrapper
    # dispatches on OpaqueParams, OpaqueVoid unpacks back to the concrete `P`, and
    # the user's `f` runs fully specialized. Lock that in so a refactor that makes
    # the unpack dynamic (e.g. losing the `P` type parameter, or a Union-typed
    # pack) is caught here rather than only showing up as a runtime slowdown.
    @testset "opaque wrapped-f call is type-stable and allocation-free" begin
        cp = concretize(deprob([1.0], LinearP(0.5)))
        u = [1.0]
        # function barrier: keeps @allocated / @inferred from capturing
        # non-const global access, so they measure the wrapped call itself.
        call(f, du, u, p, t) = f(du, u, p, t)
        measalloc(f, du, u, p, t) = (call(f, du, u, p, t); @allocated call(f, du, u, p, t))

        @test measalloc(cp.f.f, [0.0], u, cp.p, 0.0) == 0
        @test (@inferred call(cp.f.f, [0.0], u, cp.p, 0.0)) === nothing
        # the unpack resolves to the concrete payload type, not a Union
        @test @inferred(RespecializeParams.unsafe_unpack(cp.p, LinearP)) === LinearP(0.5)
        @test Base.return_types(RespecializeParams.pack_auto, (LinearP,)) ==
            [RespecializeParams.OpaqueParams]

        # AutoDePSpecialize adds no inference instability over AutoSpecialize:
        # get_concrete_problem yields the same number of (concrete) return types.
        for spec in (SciMLBase.AutoSpecialize, DiffEqBase.AutoDePSpecialize)
            prob = ODEProblem{true, spec}(linear_rhs!, [1.0], (0.0, 1.0), LinearP(0.5))
            rts = Base.return_types(DiffEqBase.get_concrete_problem, (typeof(prob), Bool))
            @test length(rts) == 1
            @test all(isconcretetype, rts)
        end
    end
end
