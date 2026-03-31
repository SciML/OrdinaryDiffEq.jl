using OrdinaryDiffEq, Test
function lorenz(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 1.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob, Tsit5(), save_idxs = 1)
@inferred solve(prob, Tsit5())
@inferred solve(prob, Tsit5(), save_idxs = 1)
@test_broken @inferred(remake(prob, u0 = Float32[1.0; 0.0; 0.0])) ==
    remake(prob, u0 = Float32[1.0; 0.0; 0.0])
@test_broken @inferred(solve(prob, Tsit5(), u0 = Float32[1.0; 0.0; 0.0])) ==
    solve(prob, Tsit5(), u0 = Float32[1.0; 0.0; 0.0])

prob = ODEProblem{true, SciMLBase.FullSpecialize}(lorenz, u0, tspan)

@inferred SciMLBase.wrapfun_iip(prob.f)
@inferred remake(prob, u0 = [1.0; 0.0; 0.0])
@inferred remake(prob, u0 = Float32[1.0; 0.0; 0.0])
@test_broken @inferred(solve(prob, Tsit5(), u0 = Float32[1.0; 0.0; 0.0])) ==
    solve(prob, Tsit5(), u0 = Float32[1.0; 0.0; 0.0])

prob = ODEProblem(lorenz, Float32[1.0; 0.0; 0.0], tspan)
@inferred solve(prob, Tsit5(), save_idxs = 1)
@test_broken @inferred(solve(prob, Tsit5(), u0 = [1.0; 0.0; 0.0])) ==
    solve(prob, Tsit5(), u0 = [1.0; 0.0; 0.0])
remake(prob, u0 = [1.0; 0.0; 0.0])

@inferred SciMLBase.wrapfun_iip(prob.f)
@test_broken @inferred(
    ODEFunction{
        isinplace(prob), SciMLBase.FunctionWrapperSpecialize,
    }(prob.f)
) ==
    ODEFunction{isinplace(prob), SciMLBase.FunctionWrapperSpecialize}(prob.f)
@inferred remake(prob, u0 = [1.0; 0.0; 0.0])
@test_broken @inferred(solve(prob, Tsit5(), u0 = [1.0; 0.0; 0.0])) ==
    solve(prob, Tsit5(), u0 = [1.0; 0.0; 0.0])

function f(du, u, p, t)
    du[1] = p.a
    return du[2] = p.b
end

const alg = Tsit5()

function solve_ode(f::F, p::P, ensemblealg; kwargs...) where {F, P}
    tspan = (0.0, 1.0)
    Δt = tspan[2] - tspan[1]
    dt = 1 / 252
    nodes = Int(ceil(Δt / dt) + 1)
    t = T = [tspan[1] + (i - 1) * dt for i in 1:nodes]

    # if I do not set {true}, prob type Any...
    prob = ODEProblem{true}(f, [0.0, 0.0], tspan, p)
    # prob = ODEProblem(f, [0., 0.], tspan, p)

    prob_func = (prob, i, repeat) -> begin
        remake(prob, tspan = (T[i + 1], t[1]))
    end

    # ensemble problem
    odes = EnsembleProblem(prob, prob_func = prob_func)

    sol = OrdinaryDiffEq.solve(
        odes, OrdinaryDiffEq.Tsit5(), ensemblealg,
        trajectories = nodes - 1, saveat = -dt;
        kwargs...
    )

    return sol
end
@inferred solve_ode(f, (a = 1, b = 1), EnsembleSerial())
@inferred solve_ode(f, (a = 1, b = 1), EnsembleThreads())
@test_broken @inferred(solve_ode(f, (a = 1, b = 1), EnsembleDistributed())) ==
    solve_ode(f, (a = 1, b = 1), EnsembleDistributed())
@test_broken @inferred(solve_ode(f, (a = 1, b = 1), EnsembleSplitThreads())) ==
    solve_ode(f, (a = 1, b = 1), EnsembleSplitThreads())
@inferred solve_ode(f, (a = 1, b = 1), EnsembleSerial(), save_idxs = 1)
@inferred solve_ode(f, (a = 1, b = 1), EnsembleThreads(), save_idxs = 1)
@test_broken @inferred(
    solve_ode(
        f, (a = 1, b = 1), EnsembleDistributed(), save_idxs = 1
    )
) ==
    solve_ode(f, (a = 1, b = 1), EnsembleDistributed(), save_idxs = 1)
@test_broken @inferred(
    solve_ode(
        f, (a = 1, b = 1), EnsembleSplitThreads(), save_idxs = 1
    )
) ==
    solve_ode(f, (a = 1, b = 1), EnsembleSplitThreads(), save_idxs = 1)

using StochasticDiffEq, Test
u0 = 1 / 2
ff(u, p, t) = u
gg(u, p, t) = u
dt = 1 // 2^(4)
tspan = (0.0, 1.0)
prob = SDEProblem(ff, gg, u0, (0.0, 1.0))
sol = solve(prob, EM(), dt = dt)
@inferred solve(prob, EM(), dt = dt)
@inferred solve(prob, EM(), dt = dt, save_idxs = 1)
