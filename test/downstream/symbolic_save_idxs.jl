using ModelingToolkit, OrdinaryDiffEq, Test

@testset "Symbolic save_idxs - Save Observed" begin
    @parameters t
    @variables a(t) b(t) c(t) d(t) e(t)

    D = Differential(t)

    eqs = [D(a) ~ a,
        D(b) ~ b,
        D(c) ~ c,
        D(d) ~ d,
        e ~ d]

    @named sys = ODESystem(eqs, t, [a, b, c, d, e], [];
                           defaults = Dict([a => 1.0,
                                               b => 1.0,
                                               c => 1.0,
                                               d => 1.0,
                                               e => 1.0]))
    sys = structural_simplify(sys)
    prob = ODEProblem(sys, [], (0, 1.0))
    prob_sym = ODEProblem(sys, [], (0, 1.0), save_idxs = [a, c, e])

    sol = solve(prob, Tsit5())
    sol_sym = solve(prob_sym, Tsit5())

    @test sol_sym[a] ≈ sol[a]
    @test sol_sym[c] ≈ sol[c]
    @test sol_sym[d] ≈ sol[d]
    @test sol_sym[e] ≈ sol[e]

    @test sol.u != sol_sym.u

    @test_throws Exception sol_sym[b]
end

@testset "Symbolic save_idxs - No observed" begin
    @parameters t
    @variables a(t) b(t) c(t) d(t) e(t)

    D = Differential(t)

    eqs = [D(a) ~ a,
        D(b) ~ b,
        D(c) ~ c,
        D(d) ~ d,
        e ~ d]

    @named sys = ODESystem(eqs, t, [a, b, c, d, e], [];
                           defaults = Dict([a => 1.0,
                                               b => 1.0,
                                               c => 1.0,
                                               d => 1.0,
                                               e => 1.0]))
    sys = structural_simplify(sys)
    prob = ODEProblem(sys, [], (0, 1.0))
    prob_sym = ODEProblem(sys, [], (0, 1.0), save_idxs = [a, c, b])

    sol = solve(prob, Tsit5())
    sol_sym = solve(prob_sym, Tsit5())

    @test sol_sym[a] ≈ sol[a]
    @test sol_sym[b] ≈ sol[b]
    @test sol_sym[c] ≈ sol[c]

    @test sol.u != sol_sym.u

    @test_throws Exception sol_sym[d]
    @test_throws Exception sol_sym[e]
end
