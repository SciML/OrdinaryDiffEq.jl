using OrdinaryDiffEqStabilizedIRK, Test, LinearAlgebra, Random
using OrdinaryDiffEqStabilizedIRK: maxeig!

@testset "Power Iteration of Runge-Kutta-Chebyshev Tests" begin
    Random.seed!(123)
    eigen_est = (integrator) -> integrator.eigen_est = 1.5e2
    for iip in [true, false], Alg in [IRKC]
        alg = Alg()
        println(typeof(alg))
        A = randn(20, 20)
        B = randn(20, 20)
        test_f1 = !iip ? (u, p, t) -> A * u : (du, u, p, t) -> mul!(du, A, u)
        test_f2 = !iip ? (u, p, t) -> B * u : (du, u, p, t) -> mul!(du, B, u)
        ff_split = SplitFunction{iip}(test_f1, test_f2)
        prob = SplitODEProblem{iip}(ff_split, randn(20, 1), (0.0, 1.0))
        integrator = init(prob, alg)
        eigm = maximum(abs.(eigvals(A)))
        maxeig!(integrator, integrator.cache)
        eigest = integrator.eigen_est
        @test eigest ≈ eigm rtol = 0.1eigm

        A = A - 1.0e2I
        test_f1 = !iip ? (u, p, t) -> A * u : (du, u, p, t) -> mul!(du, A, u)
        prob = SplitODEProblem{iip}(
            SplitFunction{iip}(test_f1, test_f2), ones(20),
            (0.0, 1.0)
        )
        @test_nowarn solve(prob, alg)
    end
end

# `f1ⱼ₋₂ = du₁` in the in-place `perform_step!` rebound the name onto `du₁`, so the
# stage shift `@.. f1ⱼ₋₂ = f1ⱼ₋₁` wrote through to `du₁` -- which has to hold
# `f1(uprev)` for the whole step. That cost the in-place method an order. The
# out-of-place method does the same assignment safely, because there it is a value
# copy, so the two paths disagreeing is the sharpest signal.
@testset "IRKC in-place matches out-of-place" begin
    A1 = [-100.0 0.0; 0.0 -50.0]
    A2 = [0.0 1.0; -1.0 0.0]
    u0 = [1.0, 1.0]
    tspan = (0.0, 0.1)
    exact = exp((A1 + A2) * 0.1) * u0

    f1! = (du, u, p, t) -> (mul!(du, A1, u); nothing)
    f2! = (du, u, p, t) -> (mul!(du, A2, u); nothing)
    f1 = (u, p, t) -> A1 * u
    f2 = (u, p, t) -> A2 * u
    prob_iip = SplitODEProblem(f1!, f2!, u0, tspan)
    prob_oop = SplitODEProblem(f1, f2, u0, tspan)

    dts = [0.1 / 2^i for i in 3:7]
    err(prob) = [
        maximum(abs.(solve(prob, IRKC(); dt = dt, adaptive = false).u[end] .- exact))
            for dt in dts
    ]
    eiip = err(prob_iip)
    eoop = err(prob_oop)

    # Same problem, same method: the two paths must agree closely.
    @test maximum(abs.(eiip .- eoop)) < 1.0e-10

    # ...and the in-place path must show the method's order, not one less.
    orders = [log2(eiip[i] / eiip[i + 1]) for i in 1:(length(dts) - 1)]
    @test minimum(orders) > 1.5
end
