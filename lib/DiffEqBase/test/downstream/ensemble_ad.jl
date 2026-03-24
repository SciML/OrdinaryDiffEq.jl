using OrdinaryDiffEq, Test
using SciMLSensitivity
using Random
using DifferentiationInterface
using ADTypes: AutoForwardDiff, AutoMooncake

# Load backends for all versions (required for DifferentiationInterface extensions)
using ForwardDiff
using Mooncake

# Version-dependent imports
if VERSION <= v"1.11"
    using Zygote
    using ADTypes: AutoZygote
end
if VERSION <= v"1.11"
    using Enzyme
    using ADTypes: AutoEnzyme
end

# Define backends based on Julia version
# ForwardDiff: All versions
# Mooncake: All versions
# Zygote: Julia <= 1.11
# Enzyme: Julia <= 1.11
function get_test_backends()
    backends = Pair{String, Any}[]
    # ForwardDiff on all versions
    push!(backends, "ForwardDiff" => AutoForwardDiff())
    # Mooncake on all versions
    push!(backends, "Mooncake" => AutoMooncake(; config = nothing))
    # Zygote only on Julia <= 1.11
    if VERSION <= v"1.11"
        push!(backends, "Zygote" => AutoZygote())
    end
    # Enzyme only on Julia <= 1.11
    if VERSION <= v"1.11"
        push!(backends, "Enzyme" => AutoEnzyme())
    end
    return backends
end

function dt!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α * x - β * x * y
    return du[2] = dy = -δ * y + γ * x * y
end

n_par = 3
Random.seed!(2)
u0 = rand(2, n_par)
u0[:, 1] = [1.0, 1.0]
tspan = (0.0, 10.0)
p = [2.2, 1.0, 2.0, 0.4]
prob_ode = ODEProblem(dt!, u0[:, 1], tspan)

function test_loss(p1, prob)
    function prob_func(prob, i, repeat)
        @show i
        return remake(prob, u0 = u0[:, i])
    end

    #define ensemble problem
    ensembleprob = EnsembleProblem(prob, prob_func = prob_func)

    u = Array(
        solve(
            ensembleprob, Tsit5(), EnsembleThreads(), trajectories = n_par,
            p = p1,
            sensealg = ForwardDiffSensitivity(),
            saveat = 0.1, dt = 0.001
        )
    )[:, end, :]
    loss = sum(u)
    return loss
end

test_loss(p, prob_ode)

# Test gradient computation with DifferentiationInterface for each backend
# Note: Mooncake and Enzyme are excluded for ensemble tests
# Mooncake: StackOverflowError in rule compilation
# Enzyme: Issues with ensemble problem differentiation
backends = get_test_backends()
backends_for_ensemble = filter(b -> b[1] ∉ ("Mooncake", "Enzyme"), backends)
for (name, backend) in backends_for_ensemble
    @testset "Ensemble test_loss gradient with $name" begin
        @time gs = DifferentiationInterface.gradient(p -> test_loss(p, prob_ode), backend, p)
        @test gs isa Vector
    end
end

### https://github.com/SciML/DiffEqFlux.jl/issues/595

function fiip(du, u, p, t)
    du[1] = dx = p[1] * u[1] - p[2] * u[1] * u[2]
    return du[2] = dy = -p[3] * u[2] + p[4] * u[1] * u[2]
end

p = [1.5, 1.0, 3.0, 1.0];
u0 = [1.0; 1.0];
prob = ODEProblem(fiip, u0, (0.0, 10.0), p)
sol = solve(prob, Tsit5())

function sum_of_solution(x)
    _prob = remake(prob, u0 = x[1:2], p = x[3:end])
    return sum(solve(_prob, Tsit5(), saveat = 0.1))
end

# Test sum_of_solution gradient with all backends
for (name, backend) in backends
    @testset "sum_of_solution gradient with $name" begin
        gs = DifferentiationInterface.gradient(sum_of_solution, backend, [u0; p])
        @test gs isa Vector
    end
end

# Testing ensemble problem with ForwardDiff
N = 3
eu0 = rand(N, 2)
ep = rand(N, 4)

ensemble_prob = EnsembleProblem(
    prob,
    prob_func = (
        prob, i, repeat,
    ) -> remake(
        prob,
        u0 = eu0[i, :],
        p = ep[i, :],
        saveat = 0.1
    )
)
esol = solve(ensemble_prob, Tsit5(), trajectories = N)

cache = Ref{Any}()

function sum_of_e_solution(p)
    ensemble_prob = EnsembleProblem(
        prob,
        prob_func = (
            prob, i, repeat,
        ) -> remake(
            prob,
            u0 = eu0[i, :],
            p = p[i, :],
            saveat = 0.1
        )
    )
    sol = solve(
        ensemble_prob, Tsit5(), EnsembleSerial(), trajectories = N, abstol = 1.0e-12,
        reltol = 1.0e-12
    )
    z = Array(sol.u[1])
    cache[] = sol.u[1].t
    return sum(z) # just test for the first solutions, gradients should be zero for others
end

sum_of_e_solution(ep)

# Test ensemble AD with multiple backends and compare results
# Note: Mooncake and Enzyme are excluded for ensemble tests (same as above)
@testset "Ensemble AD comparison across backends" begin
    # Use ForwardDiff as reference
    x_ref = DifferentiationInterface.gradient(sum_of_e_solution, AutoForwardDiff(), ep)

    for (name, backend) in backends_for_ensemble
        @testset "sum_of_e_solution gradient with $name" begin
            x = DifferentiationInterface.gradient(sum_of_e_solution, backend, ep)
            @test x ≈ x_ref
        end
    end
    @test cache[] == 0:0.1:10.0 # test prob.kwargs is forwarded
end
