# ncondition tests
using OrdinaryDiffEq, Test

function f(u, p, t)
    return 5 * u
end
u0 = [1.0, 1.0]
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
x = Ref(0)

condition = function (u, t, integrator)
    x[] += 1
    return 0
end
affect! = function (integrator) end
cb = ContinuousCallback(condition, affect!)
sol = solve(prob, Vern9(), callback = cb)
@test x[] == sol.stats.ncondition

condition = function (u, t, integrator)
    x[] += 1
    return t - 0.46
end
x[] = 0
cb = ContinuousCallback(condition, affect!)
sol = solve(prob, Vern9(), callback = cb)
@test x[] == sol.stats.ncondition

condition = function (u, t, integrator)
    x[] += 1
    return 1
end
x[] = 0
cb = ContinuousCallback(condition, affect!)
sol = solve(prob, Vern9(), callback = cb)
@test x[] == sol.stats.ncondition

condition = function (u, t, integrator)
    x[] += 1
    return true
end
x[] = 0
cb = DiscreteCallback(condition, affect!)
sol = solve(prob, Vern9(), callback = cb)
@test x[] == sol.stats.ncondition
