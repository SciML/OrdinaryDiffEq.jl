using OrdinaryDiffEqRosenbrock, SciMLBase, Test, LinearAlgebra

@testset "IMEXRKR_3_2 convergence" begin
    f1 = (du, u, p, t) -> (du .= -100 .* u)
    f2 = (du, u, p, t) -> (du .= sin(t))
    f = SplitFunction(f1, f2)
    prob = SplitODEProblem(f, [1.0], (0.0, 0.1))

    # Reference value from high-accuracy Rodas5P solve
    uref = 9.441483448553576e-4

    dts = [0.005, 0.0025, 0.00125, 0.000625]
    errs = [begin
        sol = solve(prob, IMEXRKR_3_2(), dt = dt, adaptive = false)
        abs(sol.u[end][1] - uref)
    end for dt in dts]

    orders = [log(errs[i - 1] / errs[i]) / log(dts[i - 1] / dts[i]) for i in 2:length(dts)]
    @test all(o -> o > 1.8, orders)
    @test all(o -> o < 2.3, orders)
end

@testset "IMEXRKR_3_2 out-of-place" begin
    f1 = (u, p, t) -> -100 .* u
    f2 = (u, p, t) -> [sin(t)]
    f = SplitFunction(f1, f2)
    prob = SplitODEProblem(f, [1.0], (0.0, 0.1))

    sol = solve(prob, IMEXRKR_3_2(), dt = 0.001, adaptive = false)
    @test abs(sol.u[end][1] - 9.441483448553576e-4) < 1e-6
end
