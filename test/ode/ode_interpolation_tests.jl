using OrdinaryDiffEq, Base.Test, DiffEqBase, DiffEqCallbacks

function test_ode(du, u, p, t)
    du[1] = p[1] - (1 - p[1])*u[1]
    return nothing
end

test_solution(t) = t <= 5 ? t : 5. * e^(-(t-5))

tspan = (0.,10.)
testtimes = linspace(tspan..., 1001)
condition(u, t, i) = t - 5
affect!(i) = i.p[1] = abs(1 - i.p[1])
cb = ContinuousCallback(condition, affect!)

algs = [Tsit5, Rosenbrock23, Rosenbrock32, Rodas5] ## Works for these
# algs = subtypes(OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm)

passed = fill(false, length(algs))


for (i, alg) in enumerate(algs)
    display("testing $alg")
    prob = ODEProblem(test_ode, [0.], tspan, [1.])
    sol = solve(prob, alg(); callback=cb)

    passed[i] = all(isapprox(sol(t)[1], test_solution(t); atol=0.05) for t in testtimes)
    @test passed[i]
end

any(.!(passed)) && warn("The following algorithms failed the continuous callback test: $(algs[.!(passed)])")
