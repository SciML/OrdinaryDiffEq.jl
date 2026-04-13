using DiffEqBase
using SciMLBase: @add_kwonly, add_kwonly
using LinearAlgebra, Test

@add_kwonly function f(a, b; c = 3, d = 4)
    (a, b, c, d)
end
@test f(1, 2) == (1, 2, 3, 4)
@test f(a = 1, b = 2) == (1, 2, 3, 4)
@test_throws ErrorException f()

@add_kwonly g(a, b; c = 3, d = 4) = (a, b, c, d)
@test g(1, 2) == (1, 2, 3, 4)
@test g(a = 1, b = 2) == (1, 2, 3, 4)

@add_kwonly h(; c = 3, d = 4) = (c, d)
@test h() == (3, 4)

@test_throws ErrorException add_kwonly(:(i(c = 3, d = 4) = (c, d)))

dprob = DiscreteProblem((u, p, t) -> 2u, 0.5, (0.0, 1.0))
@test remake(dprob) == dprob
@test remake(dprob; u0 = 1.0).u0 == 1.0

oprob = ODEProblem((u, p, t) -> 2u, 0.5, (0.0, 1.0))
@test_broken remake(oprob) == oprob # fails due to change to mutable struct due to === fallback
@test remake(oprob; u0 = 1.0).u0 == 1.0

sprob = SDEProblem((u, p, t) -> 2u, (u, p, t) -> 2u, 0.5, (0.0, 1.0))
@test remake(sprob) == sprob
@test remake(sprob; u0 = 1.0).u0 == 1.0

daeprob = DAEProblem((du, u, p, t) -> du - 2u, 0.5, 0.5, (0.0, 1.0))
@test remake(daeprob) == daeprob
@test remake(daeprob; u0 = 1.0).u0 == 1.0

ddeprob = DDEProblem((du, u, h, p, t) -> -2u, 0.5, (p, t) -> 0.0, (0.0, 1.0))
@test remake(ddeprob) == ddeprob
@test remake(daeprob; u0 = 1.0).u0 == 1.0

function f(du, u, p, t)
    du[1] = 0.2u[1]
    return du[2] = 0.4u[2]
end
u0 = ones(2)
tspan = (0, 1.0)

# Create a ODEProblem and test remake:
prob1 = SplitODEProblem(f, f, u0, tspan, Dict(), callback = nothing)
prob2 = @inferred remake(prob1; u0 = prob1.u0 .+ 1)
@test prob1.f === prob2.f
@test prob1.p === prob2.p
@test prob1.u0 .+ 1 ≈ prob2.u0
@test prob1.tspan == prob2.tspan
@test prob1.kwargs[:callback] === prob2.kwargs[:callback]
@test prob1.problem_type === prob2.problem_type

prob2 = @inferred remake(prob1; u0 = prob1.u0 .+ 1, callback = :test)
@test prob2.kwargs[:callback] == :test

# Test remake with SplitFunction:
prob1 = SplitODEProblem((u, p, t) -> u / 2, (u, p, t) -> 2u, 1.0, (0.0, 1.0))
prob2 = remake(
    prob1;  # prob1 is a ODEProblem
    f = remake(
        prob1.f;  # prob1.f is a SplitFunction
        f2 = (u, p, t) -> 3u
    )
)

# Test remake with NoiseProblem (a struct w/o isinplace type parameter):
struct DummyNoiseProcess <: SciMLBase.AbstractNoiseProcess{Int, 1, Nothing, true}
    dummy::Any
end
tspan1 = (0.0, 1.0)
tspan2 = (0.0, 2.0)
noise1 = NoiseProblem(DummyNoiseProcess(Dict()), tspan1);
noise2 = remake(noise1; tspan = tspan2);
@test noise1.noise === noise2.noise
@test noise1.tspan == tspan1
@test noise2.tspan == tspan2
@test noise1.tspan != noise2.tspan

# Test remake with TwoPointBVPFunction (manually defined):
f1 = SciMLBase.TwoPointBVPFunction((u, p, t) -> 1, ((u_a, p) -> 2, (u_b, p) -> 2))
@test_broken f2 = remake(f1; bc = ((u_a, p) -> 3, (u_b, p) -> 4))
@test_broken f1.bc() == 1
@test_broken f2.bc() == 2

# Testing remake for no recompile
u0 = [0; 2.0]
tspan = (0.0, 6.3)
prob = ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}((du, u, p, t) -> 2u, u0, tspan)

prob2 = remake(prob; u0 = [1; 2])
@test prob2.u0 == [1; 2]
@test prob2.f.f isa SciMLBase.FunctionWrappersWrappers.FunctionWrappersWrapper
prob2 = remake(prob; p = (1, 2))
@test remake(prob; p = (1, 2)).p == (1, 2)
@test prob2.f.f isa SciMLBase.FunctionWrappersWrappers.FunctionWrappersWrapper
SciMLBase.unwrapped_f(prob2.f)
