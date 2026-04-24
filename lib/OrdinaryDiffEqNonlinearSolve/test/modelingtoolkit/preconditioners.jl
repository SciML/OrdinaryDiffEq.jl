using OrdinaryDiffEq, LinearSolve, Test, IncompleteLU, SparseArrays
using OrdinaryDiffEqSDIRK: KenCarp47, TRBDF2
using OrdinaryDiffEqRosenbrock: Rosenbrock23, Rodas4, Rodas5

# Required due to https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::IncompleteLU.ILUFactorization{Tv, Ti}) where {Tv, Ti} = Tv

const N = 32
const xyd_brusselator = range(0, stop = 1, length = N)
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a

const iter = Ref(0)
function brusselator_2d_loop(du, u, p, t)
    global iter[] += 1
    A, B, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
        ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N),
            limit(j - 1, N)
        du[i, j, 1] = alpha * (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] - 4u[i, j, 1]) + B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] + brusselator_f(x, y, t)
        du[i, j, 2] = alpha * (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] - 4u[i, j, 2]) + A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end
    return nothing
end
p = (3.4, 1.0, 10.0, step(xyd_brusselator))

function init_brusselator_2d(xyd)
    N = length(xyd)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
        x = xyd[I[1]]
        y = xyd[I[2]]
        u[I, 1] = 22 * (y * (1 - y))^(3 / 2)
        u[I, 2] = 27 * (x * (1 - x))^(3 / 2)
    end
    return u
end
u0 = init_brusselator_2d(xyd_brusselator)
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop, u0, (0.0, 11.5), p)

# Build the Jacobian sparsity pattern for the 2D Brusselator directly.
# Each (i,j) grid point with 2 species couples to its 4 neighbors (diffusion)
# and to itself (reaction terms involving both species at the same point).
function brusselator_jac_sparsity(N)
    nvar = N * N * 2
    I = Int[]
    J = Int[]
    for i in 1:N, j in 1:N
        # Linear indices for species 1 and 2 at grid point (i,j)
        idx1 = (j - 1) * N + i          # species 1
        idx2 = idx1 + N * N             # species 2

        # Self-coupling (reaction): both species at (i,j) depend on each other
        for s in (idx1, idx2), d in (idx1, idx2)
            push!(I, d); push!(J, s)
        end

        # Diffusion stencil: 4 neighbors for each species
        for di in (-1, 1)
            ni = limit(i + di, N)
            nidx1 = (j - 1) * N + ni
            nidx2 = nidx1 + N * N
            push!(I, idx1); push!(J, nidx1)
            push!(I, idx2); push!(J, nidx2)
        end
        for dj in (-1, 1)
            nj = limit(j + dj, N)
            nidx1 = (nj - 1) * N + i
            nidx2 = nidx1 + N * N
            push!(I, idx1); push!(J, nidx1)
            push!(I, idx2); push!(J, nidx2)
        end
    end
    return sparse(I, J, ones(length(I)), nvar, nvar)
end
jac = brusselator_jac_sparsity(N)

prob_ode_brusselator_2d_sparse = ODEProblem(
    ODEFunction(
        brusselator_2d_loop,
        jac_prototype = float.(jac)
    ),
    u0, (0.0, 11.5), p
)

function incompletelu(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu(convert(AbstractMatrix, W), τ = 50.0)
    else
        Pl = Plprev
    end
    return Pl, nothing
end

iter[] = 0
sol1 = solve(
    prob_ode_brusselator_2d, KenCarp47(linsolve = KrylovJL_GMRES()),
    save_everystep = false
);
iter1 = iter[];
iter[] = 0;
sol2 = solve(
    prob_ode_brusselator_2d_sparse,
    KenCarp47(
        linsolve = KrylovJL_GMRES(), precs = incompletelu,
        concrete_jac = true
    ), save_everystep = false
);
iter2 = iter[];

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test iter2 < iter1

iter[] = 0;
sol1 = solve(
    prob_ode_brusselator_2d, Rosenbrock23(linsolve = KrylovJL_GMRES()),
    save_everystep = false
);
iter1 = iter[];
iter[] = 0;
sol2 = solve(
    prob_ode_brusselator_2d_sparse,
    Rosenbrock23(
        linsolve = KrylovJL_GMRES(), precs = incompletelu,
        concrete_jac = true
    ), save_everystep = false
);
iter2 = iter[];

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test iter2 < iter1

iter[] = 0;
sol1 = solve(
    prob_ode_brusselator_2d, Rodas4(linsolve = KrylovJL_GMRES()),
    save_everystep = false
);
iter1 = iter[];
iter[] = 0;
sol2 = solve(
    prob_ode_brusselator_2d_sparse,
    Rodas4(linsolve = KrylovJL_GMRES(), precs = incompletelu, concrete_jac = true),
    save_everystep = false
);
iter2 = iter[];

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test iter2 < iter1

iter[] = 0;
sol1 = solve(
    prob_ode_brusselator_2d, Rodas5(linsolve = KrylovJL_GMRES()),
    save_everystep = false
);
iter1 = iter[];
iter[] = 0;
sol2 = solve(
    prob_ode_brusselator_2d_sparse,
    Rodas5(linsolve = KrylovJL_GMRES(), precs = incompletelu, concrete_jac = true),
    save_everystep = false
);
iter2 = iter[];

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test iter2 < iter1

iter[] = 0;
sol1 = solve(
    prob_ode_brusselator_2d, TRBDF2(linsolve = KrylovJL_GMRES()),
    save_everystep = false
);
iter1 = iter[];
iter[] = 0;
sol2 = solve(
    prob_ode_brusselator_2d_sparse,
    TRBDF2(linsolve = KrylovJL_GMRES(), precs = incompletelu, concrete_jac = true),
    save_everystep = false
);
iter2 = iter[];

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test iter2 < iter1

iter[] = 0;
sol1 = solve(
    prob_ode_brusselator_2d, TRBDF2(linsolve = KrylovJL_GMRES()),
    save_everystep = false
);
iter1 = iter[];
iter[] = 0;
sol2 = solve(
    prob_ode_brusselator_2d_sparse,
    TRBDF2(
        linsolve = KrylovJL_GMRES(), precs = incompletelu,
        concrete_jac = true
    ), save_everystep = false
);
iter2 = iter[];

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test iter2 < iter1
