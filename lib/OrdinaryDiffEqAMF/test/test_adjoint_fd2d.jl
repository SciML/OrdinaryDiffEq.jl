using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqAMF
using SciMLOperators
using SciMLSensitivity
using LinearAlgebra
using Test

function setup_adjoint_fd2d(; A = 0.1, B = 0.1, N = 10, final_t = 1.0)
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

    # VJP: J(y) = J_op + diag(2y - 3y²), symmetric → J' = J
    function f_vjp!(dλ, λ, y, p, t)
        mul!(dλ, J_op, λ)
        @. dλ += (2 * y - 3 * y^2) * λ
        return nothing
    end

    u0 = [16 * (h * i) * (h * j) * (1 - h * i) * (1 - h * j) for i in 1:N for j in 1:N]
    tspan = (0.0, final_t)

    # Custom Kronecker AMF factors (never materialize N²×N²)
    gamma_op = ScalarOperator(
        1.0;
        update_func = (old, u, p, t; gamma = 1.0) -> gamma,
        accepted_kwargs = Val((:gamma,)),
    )
    custom_factors = (
        Base.kron(IdentityOperator(N) - gamma_op * A * D_op, IdentityOperator(N)),
        Base.kron(IdentityOperator(N), IdentityOperator(N) - gamma_op * B * D_op),
    )

    amf_func = build_amf_function(f!; jac = J_op, split = (Jx_op, Jy_op), amf_factors = custom_factors)

    # Add vjp field (build_amf_function doesn't support it yet)
    amf_func_with_vjp = ODEFunction{true, SciMLBase.FullSpecialize}(
        f!;
        jac_prototype = amf_func.jac_prototype,
        W_prototype = amf_func.W_prototype,
        sparsity = amf_func.sparsity,
        vjp = f_vjp!,
    )

    return (;
        amf_func_with_vjp, u0, tspan, J_op, Jx_op, Jy_op,
        custom_factors, N, final_t, f!, f_vjp!,
    )
end

@testset "Adjoint AMF fd2d (symmetric Jacobian)" begin
    setup = setup_adjoint_fd2d(; N = 10)
    (; amf_func_with_vjp, u0, tspan, J_op, custom_factors, N, final_t, f!, f_vjp!) = setup

    amf_prob = ODEProblem(amf_func_with_vjp, u0, tspan)

    # Forward solve
    forward_sol = solve(amf_prob, AMF(ROS34PW1a); dense = true)
    @test forward_sol.retcode == SciMLBase.ReturnCode.Success

    # Loss: L = 0.5||u(T)||², ∂L/∂u = u
    function dgdu_discrete(out, u, p, t, i)
        copyto!(out, u)
        return nothing
    end

    adj_sparsity = convert(AbstractMatrix, J_op)

    # Helper: fresh baseline adjoint problem (callback has mutable state)
    function make_baseline_prob()
        ap = SciMLSensitivity.ODEAdjointProblem(
            forward_sol,
            QuadratureAdjoint(autojacvec = false),
            ROS34PW1a(),
            [final_t],
            dgdu_discrete,
        )
        bf = ODEFunction{true, SciMLBase.FullSpecialize}(
            ap.f.f;
            jac_prototype = ap.f.jac_prototype,
            sparsity = adj_sparsity,
        )
        return ODEProblem(
            bf, ap.u0, ap.tspan, ap.p;
            callback = ap.kwargs[:callback],
        )
    end

    # Helper: fresh AMF adjoint problem
    function make_amf_prob()
        ap = SciMLSensitivity.ODEAdjointProblem(
            forward_sol,
            QuadratureAdjoint(autojacvec = false),
            ROS34PW1a(),
            [final_t],
            dgdu_discrete,
        )
        W_base = AMFOperator(; factors = custom_factors)
        W_op = cache_operator(W_base, zeros(N^2))
        af = ODEFunction{true, SciMLBase.FullSpecialize}(
            ap.f.f;
            jac_prototype = J_op,
            W_prototype = W_op,
            sparsity = adj_sparsity,
        )
        return remake(ap; f = af)
    end

    # Baseline adjoint solve (no AMF)
    adjoint_sol_baseline = solve(make_baseline_prob(), ROS34PW1a())
    @test adjoint_sol_baseline.retcode == SciMLBase.ReturnCode.Success
    dLdu0_baseline = adjoint_sol_baseline.u[end]

    # AMF adjoint solve
    adjoint_sol_amf = solve(make_amf_prob(), AMF(ROS34PW1a))
    @test adjoint_sol_amf.retcode == SciMLBase.ReturnCode.Success
    dLdu0_amf = adjoint_sol_amf.u[end]

    # AMF gradient should be close to baseline
    relerr_amf_vs_baseline = norm(dLdu0_amf - dLdu0_baseline) / norm(dLdu0_baseline)
    @test relerr_amf_vs_baseline < 1e-4
    @info "Adjoint AMF vs baseline" relerr = relerr_amf_vs_baseline

    # Finite-difference validation on a few components
    ref_func = ODEFunction{true, SciMLBase.FullSpecialize}(
        f!; jac_prototype = J_op, sparsity = adj_sparsity,
    )

    function fd_gradient_component(idx; ε = 1e-7)
        u_plus = copy(u0); u_plus[idx] += ε
        u_minus = copy(u0); u_minus[idx] -= ε
        sol_plus = solve(ODEProblem(ref_func, u_plus, tspan), ROS34PW1a(); abstol = 1e-10, reltol = 1e-10)
        sol_minus = solve(ODEProblem(ref_func, u_minus, tspan), ROS34PW1a(); abstol = 1e-10, reltol = 1e-10)
        L_plus = 0.5 * dot(sol_plus.u[end], sol_plus.u[end])
        L_minus = 0.5 * dot(sol_minus.u[end], sol_minus.u[end])
        return (L_plus - L_minus) / (2ε)
    end

    test_indices = [1, N^2 ÷ 2, N^2]
    for idx in test_indices
        fd_grad = fd_gradient_component(idx)
        adj_amf = dLdu0_amf[idx]
        relerr_vs_fd = abs(fd_grad - adj_amf) / max(abs(fd_grad), 1e-12)
        @test relerr_vs_fd < 0.01  # AMF + forward tolerance → ~1e-3
        @info "FD validation idx=$idx" fd = fd_grad adj_amf = adj_amf relerr = relerr_vs_fd
    end
end
