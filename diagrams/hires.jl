using OrdinaryDiffEq, ParameterizedFunctions, Plots, LSODA, DiffEqDevTools, Sundials
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)

gr() #gr(fmt=:png)

f_hires = @ode_def Hires begin
    dy1 = -1.71 * y1 + 0.43 * y2 + 8.32 * y3 + 0.0007
    dy2 = 1.71 * y1 - 8.75 * y2
    dy3 = -10.03 * y3 + 0.43 * y4 + 0.035 * y5
    dy4 = 8.32 * y2 + 1.71 * y3 - 1.12 * y4
    dy5 = -1.745 * y5 + 0.43 * y6 + 0.43 * y7
    dy6 = -280.0 * y6 * y8 + 0.69 * y4 + 1.71 * y5 -
          0.43 * y6 + 0.69 * y7
    dy7 = 280.0 * y6 * y8 - 1.81 * y7
    dy8 = -280.0 * y6 * y8 + 1.81 * y7
end

u0 = zeros(8)
u0[1] = 1
u0[8] = 0.0057
prob = ODEProblem(f_hires, u0, (0.0, 321.8122))

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
DIAGRAMS["HIRES-HighTol-Group1"] = plot(wp)

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
DIAGRAMS["HIRES-LowTol-Group2"] = plot(wp)
