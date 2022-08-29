using BenchmarkTools, OrdinaryDiffEq
BenchmarkTools.DEFAULT_PARAMETERS.gcsample = true
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 300
f(u,p,t) = u
prob = ODEProblem(f,1.0,(0.0,1.0))

SUITE = BenchmarkGroup()
SUITE["algs"] = BenchmarkGroup()

SUITE["algs"]["Tsit5"] = @benchmarkable solve(prob, Tsit5())
SUITE["algs"]["Vern6"] = @benchmarkable solve(prob, Vern6())
SUITE["algs"]["Vern7"] = @benchmarkable solve(prob, Vern7())
SUITE["algs"]["Vern9"] = @benchmarkable solve(prob, Vern9())
SUITE["algs"]["Rosenbrock23"] = @benchmarkable solve(prob, Rosenbrock23())
SUITE["algs"]["BS3"] = @benchmarkable solve(prob, BS3())
SUITE["algs"]["KenCarp4"] = @benchmarkable solve(prob, KenCarp4())
SUITE["algs"]["TRBDF2"] = @benchmarkable solve(prob, TRBDF2())
