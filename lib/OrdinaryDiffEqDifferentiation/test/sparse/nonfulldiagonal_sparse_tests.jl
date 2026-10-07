using OrdinaryDiffEq, SparseArrays, LinearSolve, LinearAlgebra, Test
using OrdinaryDiffEqSDIRK
using ComponentArrays
using SciMLOperators: MatrixOperator

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
        return MatrixOperator(1 / dx * M)
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
        return MatrixOperator(1 / dx^2 * M)
    end

    function extender(N)
        dx = 1 / (N + 1)
        diag = ones(N)
        lower = spzeros(Float64, N)
        upper = spzeros(Float64, N)
        lower[1] = 1.0
        upper[end] = 1.0
        M = vcat(
            transpose(lower),
            sparse(diagm(diag)),
            transpose(upper)
        )
        return MatrixOperator(1 / dx^2 * M)
    end

    bc_handler = extender(N)

    ∇ = first_deriv(N) * bc_handler
    Δ = second_deriv(N) * bc_handler

    bc_x = zeros(Real, N)
    bc_xx = zeros(Real, N)

    function timedifferentialclosure!(du, u, p, t)
        (;
            α, D, v, k_p, V_c, Q_l, Q_r, V_b,
            S, Lm, Dm, V_v,
        ) = p

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
        return
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
    V_v = 18.0
)

r_space = collect(range(0.0, 2.0, length = 15))
computeparams = (;
    Δr = r_space[2],
    r_space,
    countorderapprox = 2,
)
parameters = (;
    prior,
    compute = computeparams,
)

dudt = enclosethetimedifferential(parameters)
IC = ones(length(r_space) + 3)
odeprob = ODEProblem(
    dudt,
    IC,
    (0, 2.1),
    parameters.prior
);
du0 = copy(odeprob.u0);
# Hardcoded sparsity pattern for 15 spatial points + 3 state variables (18x18 matrix)
I = [1, 2, 16, 18, 1, 2, 3, 18, 2, 3, 4, 18, 3, 4, 5, 18, 4, 5, 6, 18, 5, 6, 7, 18, 6, 7, 8, 18, 7, 8, 9, 18, 8, 9, 10, 18, 9, 10, 11, 18, 10, 11, 12, 18, 11, 12, 13, 18, 12, 13, 14, 18, 13, 14, 15, 18, 14, 15, 17, 18, 1, 16, 17, 15, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18]
J = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 16, 16, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18]
jac_sparsity = sparse(I, J, ones(Bool, length(I)), 18, 18);
f = ODEFunction(
    dudt;
    jac_prototype = float.(jac_sparsity)
);
sparseodeprob = ODEProblem(
    f,
    odeprob.u0,
    (0, 2.1),
    parameters.prior
);

solve(odeprob, TRBDF2());
solve(sparseodeprob, TRBDF2());
solve(sparseodeprob, Rosenbrock23(linsolve = KLUFactorization()));
solve(sparseodeprob, KenCarp47(linsolve = KrylovJL_GMRES()));

@testset "Sparse Jacobian caches are initialized with stored zeros" begin
    function sparse_cache_f!(du, u, p, t)
        return du .= u
    end

    N = 8
    jac_prototype = sparse(1:N, 1:N, ones(N), N, N)
    prob = ODEProblem(
        ODEFunction(sparse_cache_f!; jac_prototype),
        ones(N),
        (0.0, 1.0)
    )
    integ = init(prob, Rodas5P())

    @test all(iszero, nonzeros(integ.cache.J))
    @test all(iszero, nonzeros(integ.cache.W))
end

@testset "Sparse W zero-init enables LinearSolve nonstructural-zero Auto reduction" begin
    n = 4
    jp = spdiagm(-1 => ones(n - 1), 0 => ones(n), 1 => ones(n - 1))
    f!(du, u, p, t) = (du .= (-u); nothing)
    function jac!(J, u, p, t)
        nonzeros(J) .= 0.0
        for i in 1:n
            J[i, i] = -1.0
        end
        return nothing
    end
    function sparse_reduction(lc)
        return lc.cacheval isa LinearSolve.DefaultLinearSolverInit ?
            lc.cacheval.sparse_reduction : lc.sparse_reduction
    end

    prob = ODEProblem(ODEFunction(f!; jac = jac!, jac_prototype = jp), ones(n), (0.0, 1.0))
    integ = init(prob, Rodas5P(linsolve = PureKLUFactorization()))
    red = sparse_reduction(integ.cache.linsolve)
    @test all(iszero, nonzeros(integ.cache.J))
    @test all(iszero, nonzeros(integ.cache.W))
    @test !red.active
    @test red.pending

    step!(integ)
    red = sparse_reduction(integ.cache.linsolve)
    @test count(iszero, nonzeros(integ.cache.W)) >= 1
    @test !red.pending
    @test red.active
    @test red.nstart_zeros > 0
end
