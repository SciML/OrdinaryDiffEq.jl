using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqAMF
using SciMLOperators
using LinearAlgebra
using Test

function setup_fd2d_problem(; A = 0.1, B = 0.1, C = 0.0, N = 40, final_t = 1.0)
    @assert iszero(C) "mixed derivative term not handled in this regression test"

    h = 1 / (N + 1)
    D = 1 / h^2 * Tridiagonal(ones(N - 1), -2 * ones(N), ones(N - 1))
    D_op = MatrixOperator(D)

    Jx_op = A * Base.kron(D_op, IdentityOperator(N))
    Jy_op = B * Base.kron(IdentityOperator(N), D_op)
    J_op = cache_operator(Jx_op + Jy_op, zeros(N^2))

    function g!(du, u, p, t)
        @. du += u^2 * (1 - u) + exp(t)
        return nothing
    end

    function f!(du, u, p, t)
        mul!(du, J_op, u)
        g!(du, u, p, t)
        return nothing
    end

    u0 = [16 * (h * i) * (h * j) * (1 - h * i) * (1 - h * j) for i in 1:N for j in 1:N]
    tspan = (0.0, final_t)

    ref_func = ODEFunction{true, SciMLBase.FullSpecialize}(f!; jac_prototype = J_op)
    ref_prob = ODEProblem(ref_func, u0, tspan)

    amf_func = build_amf_function(f!; jac = J_op, split = (Jx_op, Jy_op))
    amf_prob = ODEProblem(amf_func, u0, tspan)

    gamma_op = ScalarOperator(
        1.0;
        update_func = (old, u, p, t; gamma = 1.0) -> gamma,
        accepted_kwargs = Val((:gamma,)),
    )
    custom_factors = (
        Base.kron(IdentityOperator(N) - gamma_op * A * D_op, IdentityOperator(N)),
        Base.kron(IdentityOperator(N), IdentityOperator(N) - gamma_op * B * D_op),
    )
    custom_amf_func = build_amf_function(f!; jac = J_op, split = (Jx_op, Jy_op), amf_factors = custom_factors)
    custom_amf_prob = ODEProblem(custom_amf_func, u0, tspan)

    exact_w_func = build_amf_function(f!; jac = J_op)
    exact_w_prob = ODEProblem(exact_w_func, u0, tspan)

    return ref_prob, amf_prob, custom_amf_prob, exact_w_prob
end

@testset "AMF finite-difference 2D regression" begin
    ref_prob, amf_prob, custom_amf_prob, exact_w_prob = setup_fd2d_problem()

    @test !isa(exact_w_prob.f.W_prototype, MatrixOperator)

    ref_sol = solve(ref_prob, ROS34PW1a(); abstol = 1.0e-10, reltol = 1.0e-10)
    amf_sol = solve(amf_prob, AMF(ROS34PW1a); abstol = 1.0e-8, reltol = 1.0e-8)
    custom_amf_sol = solve(custom_amf_prob, AMF(ROS34PW1a); abstol = 1.0e-8, reltol = 1.0e-8)
    exact_w_sol = solve(exact_w_prob, AMF(ROS34PW1a); abstol = 1.0e-10, reltol = 1.0e-10)

    @test ref_sol.retcode == SciMLBase.ReturnCode.Success
    @test amf_sol.retcode == SciMLBase.ReturnCode.Success
    @test custom_amf_sol.retcode == SciMLBase.ReturnCode.Success
    @test exact_w_sol.retcode == SciMLBase.ReturnCode.Success

    relerr = norm(amf_sol.u[end] - ref_sol.u[end]) / norm(ref_sol.u[end])
    @test relerr < 1.0e-6

    custom_relerr = norm(custom_amf_sol.u[end] - ref_sol.u[end]) / norm(ref_sol.u[end])
    @test custom_relerr < 1.0e-6

    exact_w_relerr = norm(exact_w_sol.u[end] - ref_sol.u[end]) / norm(ref_sol.u[end])
    @test exact_w_relerr < 1.0e-10
end

@testset "AMF operator validation" begin
    D2 = MatrixOperator(Diagonal(ones(2)))
    D3 = MatrixOperator(Diagonal(ones(3)))

    @test_throws ArgumentError build_amf_function((du, u, p, t) -> nothing; jac = D2, split = ())
    @test_throws ArgumentError build_amf_function((du, u, p, t) -> nothing; jac = D2, split = (D3,))
    @test_throws ArgumentError build_amf_function((du, u, p, t) -> nothing; jac = D3, split = (D2, D3))
    @test_throws ArgumentError build_amf_function((du, u, p, t) -> nothing; jac = D2, amf_factors = ())
    @test_throws ArgumentError build_amf_function((du, u, p, t) -> nothing; jac = D2, amf_factors = (D3,))
    @test_throws ArgumentError build_amf_function((du, u, p, t) -> nothing; jac = D2, split = (D2,), amf_factors = (D2, D2))
end
