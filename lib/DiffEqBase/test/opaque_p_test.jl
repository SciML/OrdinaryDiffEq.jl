using DiffEqBase
using SciMLBase
using RespecializeParams
using Test

# An `isbits` parameter struct. `get_concrete_problem` on an AutoSpecialize
# in-place ODEFunction with such a `p` should route through the opaque path:
# the returned problem has `prob.p :: OpaqueParams` and `prob.f.f` is a
# FunctionWrappersWrapper whose signature carries `OpaqueParams` in the p
# slot.

struct LinearP
    k::Float64
end

function linear_rhs!(du, u, p::LinearP, t)
    @inbounds du[1] = -p.k * u[1]
    return nothing
end

# Trigger the same code path the high-level `solve(...)` does, without
# requiring a concrete integrator package as a test dep.
concretize(prob, alg = nothing) = DiffEqBase.get_concrete_problem(prob, true; alg = alg)

@testset "opaque-p hook" begin
    @testset "AutoSpecialize opaque-ifies isbits p" begin
        prob = ODEProblem(linear_rhs!, [1.0], (0.0, 1.0), LinearP(0.5))
        cp = concretize(prob)
        @test cp.p isa RespecializeParams.OpaqueParams
        # The unpack-then-call path works on the wrapped rhs:
        du = [0.0]
        cp.f(du, [1.0], cp.p, 0.0)
        @test du[1] ≈ -0.5
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
        prob = ODEProblem(f_null!, [1.0], (0.0, 1.0))
        cp = concretize(prob)
        @test cp.p isa SciMLBase.NullParameters
    end

    @testset "Non-isbits p skipped" begin
        f_vec!(du, u, p::Vector{Float64}, t) = (@inbounds du[1] = -p[1] * u[1]; nothing)
        prob = ODEProblem(f_vec!, [1.0], (0.0, 1.0), [0.5])
        cp = concretize(prob)
        @test cp.p isa Vector{Float64}
    end

    @testset "two LinearP problems share the wrapped-f type" begin
        prob_a = ODEProblem(linear_rhs!, [1.0], (0.0, 1.0), LinearP(0.5))
        prob_b = ODEProblem(linear_rhs!, [2.0], (0.0, 1.0), LinearP(1.5))
        cp_a = concretize(prob_a)
        cp_b = concretize(prob_b)
        @test typeof(cp_a.f) === typeof(cp_b.f)
        @test typeof(cp_a.p) === typeof(cp_b.p) === RespecializeParams.OpaqueParams
    end

    @testset "unsafe_unpack of the packed p yields the original value" begin
        p = LinearP(0.5)
        prob = ODEProblem(linear_rhs!, [1.0], (0.0, 1.0), p)
        cp = concretize(prob)
        @test RespecializeParams.unsafe_unpack(cp.p, LinearP) === p
    end
end
