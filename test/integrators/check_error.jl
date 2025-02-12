using OrdinaryDiffEq, Test

f_ec(u, p, t) = exp(u)
u0 = 0.0 # explosion time is 1.0
tspan = (0.0, 10.0)
prob = ODEProblem(f_ec, u0, tspan)
options = [:reltol => 1e-8, :abstol => 1e-8, :verbose => false]
desired_codes = (ReturnCode.MaxIters, ReturnCode.Unstable)

# Test that sol.retcode is set to the correct value by various ways to
# invoke integrator.

sol = solve(prob, Tsit5(); options...)
@test sol.retcode in desired_codes

integrator = init(prob, Tsit5(); options...)
solve!(integrator)
@test integrator.sol.retcode in desired_codes

integrator = init(prob, Tsit5(); options...)
for _ in integrator
end
@test integrator.sol.retcode in desired_codes

integrator = init(prob, Tsit5(); options...)
step!(integrator, 10.0)
@test integrator.sol.retcode in desired_codes

# Test check_error
integrator = init(prob, Tsit5(); options...)
step!(integrator)
@test check_error(integrator) == ReturnCode.Success
ok = false
for i in 1:(integrator.opts.maxiters)
    step!(integrator)
    if check_error(integrator) in desired_codes
        global ok = true
        # @show i
        break
    end
end
@test ok

let
    function f!(out, u, _, t)
        out[1] = u[1] + 1 - sin(t)
    end
    mprob = ODEProblem(ODEFunction(f!, mass_matrix = [0.0;;]), [0.0], (0, 2.0))
    @test solve(mprob, Rosenbrock23()).retcode == ReturnCode.Success
end

@testset "Callbacks shouldn't disable error checking" begin
    callback = ContinuousCallback((u, t, integ) -> t - prevfloat(0.5), Returns(nothing))
    prob = ODEProblem((u, p, t) -> u, 0.0, (0.0, 1); tstops = [0.5], callback)
    sol = solve(prob, FBDF(), maxiters = 30)
    @test sol.stats.naccept + sol.stats.nreject <= 30
    @test_broken sol.retcode = ReturnCode.Success
end

@test_throws ArgumentError solve(prob, Euler(), dt = 0.1, adaptive = true)
@test_throws ArgumentError solve(prob, Euler())
