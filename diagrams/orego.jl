using OrdinaryDiffEq, ParameterizedFunctions, Plots, DiffEqDevTools, Sundials, LSODA
gr() #gr(fmt=:png)
f = @ode_def Orego begin
    dy1 = p1 * (y2 + y1 * (1 - p2 * y1 - y2))
    dy2 = (y3 - (1 + y1) * y2) / p1
    dy3 = p3 * (y1 - y3)
end p1 p2 p3

p = [77.27, 8.375e-6, 0.161]
prob = ODEProblem(f, [1.0, 2.0, 3.0], (0.0, 30.0), p)
sol = solve(prob, Rodas5(), abstol = 1 / 10^14, reltol = 1 / 10^14)
test_sol = TestSolution(sol)

# group 1
setups = [Dict(:alg => Rosenbrock23()),
    Dict(:alg => TRBDF2()),
    Dict(:alg => CVODE_BDF()),
    Dict(:alg => lsoda())]

# High Tolerance
abstols = 1 ./ 10 .^ (5:8)
reltols = 1 ./ 10 .^ (1:4);

wp = WorkPrecisionSet(prob, abstols, reltols, setups; dense = false, verbose = false,
                      appxsol = test_sol, maxiters = Int(1e5), error_estimate = :l2)
DIAGRAMS["Orego-HighTol-Group1"] = plot(wp)

# group 2
setups = [Dict(:alg => Rodas4()),
    Dict(:alg => Rodas5()),
    Dict(:alg => Rodas5P()),
    Dict(:alg => KenCarp4()),
    Dict(:alg => RadauIIA5()),
    Dict(:alg => CVODE_BDF()),
    Dict(:alg => lsoda())]

# Low Tolerance
abstols = 1 ./ 10 .^ (7:13)
reltols = 1 ./ 10 .^ (4:10)

wp = WorkPrecisionSet(prob, abstols, reltols, setups; verbose = false,
                      dense = false, appxsol = test_sol, maxiters = Int(1e5),
                      error_estimate = :l2)
DIAGRAMS["Orego-LowTol-Group2"] = plot(wp)
