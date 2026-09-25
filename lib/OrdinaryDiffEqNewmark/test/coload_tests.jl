# Loading a sublibrary that carries AD preparation (SDIRK, Rosenbrock, or the
# umbrella) routes every solve through prepare_alg -> remake, which rebuilds
# the algorithm from its field names as keyword arguments. GeneralizedAlpha
# only accepted its parameters positionally, so any co-loaded solve threw a
# MethodError before the keyword constructor existed.
using OrdinaryDiffEqNewmark, OrdinaryDiffEqSDIRK, Test, RecursiveArrayTools
using SciMLBase: ReturnCode

f1_iip!(dv, v, u, p, t) = (dv .= -u)
f2_iip!(du, v, u, p, t) = (du .= v)
f1_oop(v, u, p, t) = -u
f2_oop(v, u, p, t) = v

@testset "Solves succeed with OrdinaryDiffEqSDIRK co-loaded" begin
    prob_iip = DynamicalODEProblem(
        DynamicalODEFunction(f1_iip!, f2_iip!), ones(2), zeros(2), (0.0, 1.0)
    )
    prob_oop = DynamicalODEProblem(
        DynamicalODEFunction(f1_oop, f2_oop), ones(2), zeros(2), (0.0, 1.0)
    )
    algs = (
        NewmarkBeta(),
        GeneralizedAlpha(rho_inf = 0.9),
        GeneralizedAlpha(alpha_hht = -0.1),
        GeneralizedAlpha(0.1, 0.2, 0.25, 0.6),
    )
    for prob in (prob_iip, prob_oop), alg in algs
        sol = solve(prob, alg, dt = 0.1)
        @test sol.retcode == ReturnCode.Success
    end
end
