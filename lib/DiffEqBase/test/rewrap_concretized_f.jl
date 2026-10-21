using DiffEqBase, SciMLBase, ForwardDiff, Test
import FunctionWrappersWrappers

struct ChunkedForwardDiffAlgorithm{CS} <: SciMLBase.AbstractODEAlgorithm end
SciMLBase.forwarddiffs_model(::ChunkedForwardDiffAlgorithm) = true
SciMLBase.forwarddiff_chunksize(::ChunkedForwardDiffAlgorithm{CS}) where {CS} = Val(CS)

struct NonForwardDiffAlgorithm <: SciMLBase.AbstractODEAlgorithm end

function decay!(du, u, p, t)
    du .= -p[1] .* u
    return nothing
end

concretize(prob, alg) = DiffEqBase.get_concrete_problem(prob, true; alg)

const DualTag = ForwardDiff.Tag{DiffEqBase.OrdinaryDiffEqTag, Float64}

function has_jacobian_signature(w, ::Val{CS}) where {CS}
    D = Vector{ForwardDiff.Dual{DualTag, Float64, CS}}
    return any(fw -> fieldtype(typeof(fw).parameters[2], 2) === D, w.fw)
end

# Calls `f` with the chunk-`CS` duals of a ForwardDiff Jacobian and returns the result.
function call_with_duals(prob, ::Val{CS}) where {CS}
    seeds = ntuple(i -> ntuple(j -> Float64(i == j), CS), CS)
    u = [ForwardDiff.Dual{DualTag}(prob.u0[i], seeds[i]...) for i in 1:CS]
    du = similar(u)
    prob.f(du, u, prob.p, 0.0)
    return du
end

@testset "a concretized f is rewrapped when its signatures do not cover the algorithm" begin
    for spec in (SciMLBase.AutoSpecialize, SciMLBase.AutoDespecialize)
        prob = ODEProblem{true, spec}(decay!, [1.0, 2.0, 3.0], (0.0, 1.0), [0.5])
        fresh3 = concretize(prob, ChunkedForwardDiffAlgorithm{3}())

        for first_alg in (ChunkedForwardDiffAlgorithm{1}(), NonForwardDiffAlgorithm())
            @testset "$(nameof(spec)) after $(nameof(typeof(first_alg)))" begin
                first = concretize(prob, first_alg)
                @test first.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper

                again = concretize(first, ChunkedForwardDiffAlgorithm{3}())
                @test typeof(again.f) === typeof(fresh3.f)
                @test has_jacobian_signature(again.f.f, Val(3))
                du = call_with_duals(again, Val(3))
                @test ForwardDiff.value.(du) == [-0.5, -1.0, -1.5]
                @test ForwardDiff.partials.(du) == [
                    ForwardDiff.Partials((-0.5, 0.0, 0.0)),
                    ForwardDiff.Partials((0.0, -0.5, 0.0)),
                    ForwardDiff.Partials((0.0, 0.0, -0.5)),
                ]
                @test again.f.f.cache_storage.cached === nothing
            end
        end

        # A wrapper that covers the algorithm is reused as it is.
        @test concretize(fresh3, ChunkedForwardDiffAlgorithm{3}()).f.f === fresh3.f.f
        @test concretize(fresh3, NonForwardDiffAlgorithm()).f.f === fresh3.f.f
    end
end
