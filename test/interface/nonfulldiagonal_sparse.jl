using OrdinaryDiffEq, SparseArrays, LinearSolve, LinearAlgebra
using ComponentArrays

function enclosethetimedifferential(parameters::NamedTuple)::Function
    @info "Enclosing the time differential"

    (; Δr, r_space, countorderapprox) = parameters.compute
    N = length(r_space)

    function first_deriv(N)
        dx = 1 / (N + 1)
        du = -1 * ones(N - 1) # off diagonal
        du2 = ones(N - 1) # off diagonal
        diag = zeros(N)
        lower = spzeros(Float64, N)
        upper = spzeros(Float64, N)
        lower[1] = -1.0
        upper[end] = 1.0
        M = hcat(lower, sparse(diagm(-1 => du, 0 => diag, 1 => du2)), upper)
        MatrixOperator(1 / dx * M)
    end

    function second_deriv(N)
        dx = 1 / (N + 1)
        du = ones(N - 1) # off diagonal
        du2 = ones(N - 1) # off diagonal
        diag = -2 * ones(N)
        lower = spzeros(Float64, N)
        upper = spzeros(Float64, N)
        lower[1] = 1.0
        upper[end] = 1.0
        M = hcat(lower, sparse(diagm(-1 => du, 0 => diag, 1 => du2)), upper)
        MatrixOperator(1 / dx^2 * M)
    end

    function extender(N)
        dx = 1 / (N + 1)
        diag = ones(N)
        lower = spzeros(Float64, N)
        upper = spzeros(Float64, N)
        lower[1] = 1.0
        upper[end] = 1.0
        M = vcat(transpose(lower),
            sparse(diagm(diag)),
            transpose(upper))
        MatrixOperator(1 / dx^2 * M)
    end

    bc_handler = extender(N)

    ∇ = first_deriv(N) * bc_handler
    Δ = second_deriv(N) * bc_handler

    bc_x = zeros(Real, N)
    bc_xx = zeros(Real, N)

    function timedifferentialclosure!(du, u, p, t)
        (; α, D, v, k_p, V_c, Q_l, Q_r, V_b,
        S, Lm, Dm, V_v) = p

        c = u[1:(end - 3)]
        c_v = u[end - 2]
        c_c = u[end - 1]
        c_b = u[end]

        J_B0 = (Dm / Lm) * (α * c_v - c[1])
        J_BL = (Dm / Lm) * (c[end] - α * c_c)
        grad_0 = (v ./ D) .* c[1] .- J_B0 ./ D
        grad_L = (v ./ D) .* c[end] .- J_BL ./ D

        bc_x[1] = grad_0 / 2
        bc_x[end] = grad_L / 2
        grad_c = ∇ * c + bc_x

        bc_xx[1] = -grad_0 / Δr
        bc_xx[end] = grad_L / Δr
        Lap_c = Δ * c + bc_xx

        C = sum(Δr .* S * (k_p * (c .- c_b)))

        dc_dt = D * Lap_c - v * grad_c .- k_p * (c .- c_b)
        du[1:(end - 3)] = dc_dt[1:end]

        dcv_dt = -S * J_B0 / V_v - (Q_l / V_v) * c_v
        du[end - 2] = dcv_dt

        dcc_dt = S * α * J_BL / V_c + (Q_l / V_c) * c_v - (Q_l / V_c) * c_c
        du[end - 1] = dcc_dt

        dcb_dt = (Q_l / V_b) * c_c + C / V_b
        du[end] = dcb_dt
    end

    return timedifferentialclosure!
end

prior = ComponentArray(;
    α = 0.2,
    D = 0.46,
    v = 0.0,
    k_p = 0.0,
    V_c = 18,
    Q_l = 20,
    Q_r = 3.6,
    V_b = 1490,
    S = 52,
    Lm = 0.05,
    Dm = 0.046,
    V_v = 18.0)

r_space = collect(range(0.0, 2.0, length = 15))
computeparams = (Δr = r_space[2],
    r_space = r_space,
    countorderapprox = 2)
parameters = (prior = prior,
    compute = computeparams)

dudt = enclosethetimedifferential(parameters)
IC = ones(length(r_space) + 3)
odeprob = ODEProblem(dudt,
    IC,
    (0, 2.1),
    parameters.prior);
du0 = copy(odeprob.u0);
# Hardcoded sparsity pattern for 15 spatial points + 3 state variables (18x18 matrix)
# Previously computed using: Symbolics.jacobian_sparsity((du, u) -> dudt(du, u, parameters.prior, 0.0), du0, odeprob.u0)
# This avoids the dependency on Symbolics in tests
I = [1, 2, 16, 18, 1, 2, 3, 18, 2, 3, 4, 18, 3, 4, 5, 18, 4, 5, 6, 18, 5, 6, 7, 18, 6, 7, 8, 18, 7, 8, 9, 18, 8, 9, 10, 18, 9, 10, 11, 18, 10, 11, 12, 18, 11, 12, 13, 18, 12, 13, 14, 18, 13, 14, 15, 18, 14, 15, 17, 18, 1, 16, 17, 15, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18]
J = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18]
jac_sparsity = sparse(I, J, ones(Bool, length(I)), 18, 18);
f = ODEFunction(dudt;
    jac_prototype = float.(jac_sparsity));
sparseodeprob = ODEProblem(f,
    odeprob.u0,
    (0, 2.1),
    parameters.prior);

solve(odeprob, TRBDF2());
solve(sparseodeprob, TRBDF2());
solve(sparseodeprob, Rosenbrock23(linsolve = KLUFactorization()));
solve(sparseodeprob, KenCarp47(linsolve = KrylovJL_GMRES()));
