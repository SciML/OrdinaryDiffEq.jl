using OrdinaryDiffEq, Test

f_ec(u,p,t) = exp(u)
u0 = 0.0 # explosion time is 1.0
tspan = (0.0, 10.0)
prob = ODEProblem(f_ec,u0,tspan)
options = [:reltol => 1e-8, :abstol => 1e-8, :verbose => false]
desired_code = :DtLessThanMin

# Test that sol.retcode is set to the correct value by various ways to
# invoke integrator.

sol = solve(prob,Tsit5(); options...)
@test sol.retcode == desired_code

integrator = init(prob,Tsit5(); options...)
solve!(integrator)
@test integrator.sol.retcode == desired_code

integrator = init(prob,Tsit5(); options...)
for _ in integrator end
@test integrator.sol.retcode == desired_code

integrator = init(prob,Tsit5(); options...)
step!(integrator, 10.0)
@test integrator.sol.retcode == desired_code

# Test check_error
integrator = init(prob,Tsit5(); options...)
step!(integrator)
@test check_error(integrator) == :Success
ok = false
for i in 1:integrator.opts.maxiters
  step!(integrator)
  if check_error(integrator) == desired_code
    global ok = true
    # @show i
    break
  end
end
@test ok
