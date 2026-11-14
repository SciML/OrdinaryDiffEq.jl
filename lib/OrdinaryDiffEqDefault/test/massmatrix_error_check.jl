using OrdinaryDiffEqDefault, SciMLBase, LinearAlgebra, Test

# A mass-matrix problem starts the default algorithm on its stiff branch and is never
# switched. The integrator's error checks (`maxiters`, `dtmin`, NaN `dt`, `unstable_check`)
# must stay on there: with a large Jacobian norm the stiffness test reports "stiff", which
# used to turn `do_error_check` off because the switcher still believed it was on the
# nonstiff branch.
function dae!(du, u, p, t)
    k, c = p
    du[1] = c - u[1]
    du[2] = k * (u[2]^3 + u[2] - u[1])
    return nothing
end
fun = ODEFunction(dae!; mass_matrix = Diagonal([1.0, 0.0]))
u0 = [1.0, 0.6823278038280184]  # consistent: y^3 + y = 1

@testset "maxiters on a mass-matrix problem, k = $k" for k in (1.0, 1.0e5)
    prob = ODEProblem(fun, u0, (0.0, 10.0), [k, 2.0])
    for alg in (DefaultODEAlgorithm(), DefaultODEAlgorithm(stiffalgfirst = true))
        sol = solve(prob, alg; maxiters = 3)
        @test sol.retcode == ReturnCode.MaxIters
        @test length(sol.t) == 4
    end
    # without the cap the solve succeeds and never leaves the stiff branch
    sol = solve(prob, DefaultODEAlgorithm())
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[end][1] ≈ 2.0 rtol = 1.0e-3
    @test allequal(sol.alg_choice)

    # `reinit!` resets the switcher, and the checks must stay on after it
    integ = init(prob, DefaultODEAlgorithm(); maxiters = 3)
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.MaxIters
    reinit!(integ)
    solve!(integ)
    @test integ.sol.retcode == ReturnCode.MaxIters
end
