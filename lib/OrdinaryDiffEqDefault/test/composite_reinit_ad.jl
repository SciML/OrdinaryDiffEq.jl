using OrdinaryDiffEqDefault, ForwardDiff, Test

# u(1) = exp(-p), so du(1)/dp = -exp(-1) at p = 1
function reused_decay(p)
    integ = init(
        ODEProblem((u, p, t) -> -p * u, 1.0, (0.0, 1.0), p), DefaultODEAlgorithm();
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    solve!(integ)
    reinit!(integ)
    return solve!(integ).u[end]
end
@test reused_decay(1.0) ≈ exp(-1) rtol = 1.0e-6
@test ForwardDiff.derivative(reused_decay, 1.0) ≈ -exp(-1) rtol = 1.0e-6
