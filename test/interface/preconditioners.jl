using OrdinaryDiffEq, LinearSolve, Test, IncompleteLU, ModelingToolkit
import AlgebraicMultigrid

# Required due to https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::IncompleteLU.ILUFactorization{Tv, Ti}) where {Tv, Ti} = Tv
Base.eltype(::AlgebraicMultigrid.Preconditioner) = Float64

const N = 32
const xyd_brusselator = range(0, stop = 1, length = N)
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0
limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a

const iter = Ref(0)
function brusselator_2d_loop(du, u, p, t)
    global iter[] += 1

    u = reshape(u, N, N, 2)
    du = reshape(du, N, N, 2)

    A, B, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x, y = xyd_brusselator[I[1]], xyd_brusselator[I[2]]
        ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N),
                             limit(j - 1, N)
        du[i, j, 1] = alpha * (u[im1, j, 1] + u[ip1, j, 1] + u[i, jp1, 1] + u[i, jm1, 1] -
                       4u[i, j, 1]) +
                      B + u[i, j, 1]^2 * u[i, j, 2] - (A + 1) * u[i, j, 1] +
                      brusselator_f(x, y, t)
        du[i, j, 2] = alpha * (u[im1, j, 2] + u[ip1, j, 2] + u[i, jp1, 2] + u[i, jm1, 2] -
                       4u[i, j, 2]) +
                      A * u[i, j, 1] - u[i, j, 1]^2 * u[i, j, 2]
    end

    vec(du)
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
    u
end
u0 = init_brusselator_2d(xyd_brusselator) |> vec
prob_ode_brusselator_2d = ODEProblem(brusselator_2d_loop, u0, (0.0, 11.5), p)

du0 = copy(u0)
jac = ModelingToolkit.Symbolics.jacobian_sparsity((du, u) -> brusselator_2d_loop(du, u, p,
                                                                                 0.0), du0,
                                                  u0)

prob_ode_brusselator_2d_sparse = ODEProblem(ODEFunction(brusselator_2d_loop,
                                                        jac_prototype = float.(jac)),
                                            u0, (0.0, 11.5), p)

function incompletelu(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = ilu(convert(AbstractMatrix, W), τ = 50.0)
    else
        Pl = Plprev
    end
    Pl, nothing
end

function algebraicmultigrid(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(convert(AbstractMatrix,
                                                                                        W)))
    else
        Pl = Plprev
    end
    Pl, nothing
end

function algebraicmultigrid2(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        A = convert(AbstractMatrix, W)
        Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(A,
                                                                                presmoother = AlgebraicMultigrid.Jacobi(rand(size(A,
                                                                                                                                  1))),
                                                                                postsmoother = AlgebraicMultigrid.Jacobi(rand(size(A,
                                                                                                                                   1)))))
    else
        Pl = Plprev
    end
    Pl, nothing
end

iter[] = 0
sol1 = solve(prob_ode_brusselator_2d, KenCarp47(linsolve = KrylovJL_GMRES()),
             save_everystep = false);
iter1 = iter[];
iter[] = 0;
sol2 = solve(prob_ode_brusselator_2d_sparse,
             KenCarp47(linsolve = KrylovJL_GMRES(), precs = incompletelu,
                       concrete_jac = true), save_everystep = false);
iter2 = iter[];
iter[] = 0;
sol3 = solve(prob_ode_brusselator_2d_sparse,
             KenCarp47(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid,
                       concrete_jac = true), save_everystep = false);
iter3 = iter[];
iter[] = 0;
sol4 = solve(prob_ode_brusselator_2d_sparse,
             KenCarp47(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid2,
                       concrete_jac = true), save_everystep = false);
iter4 = iter[];
iter[] = 0;

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test sol3.retcode === ReturnCode.Success
@test sol4.retcode === ReturnCode.Success

@test iter2 < iter1
@test iter3 < iter1
@test iter4 < iter1

sol1 = solve(prob_ode_brusselator_2d, Rosenbrock23(linsolve = KrylovJL_GMRES()),
             save_everystep = false);
iter1 = iter[];
iter[] = 0;
sol2 = solve(prob_ode_brusselator_2d_sparse,
             Rosenbrock23(linsolve = KrylovJL_GMRES(), precs = incompletelu,
                          concrete_jac = true), save_everystep = false);
iter2 = iter[];
iter[] = 0;
sol3 = solve(prob_ode_brusselator_2d_sparse,
             Rosenbrock23(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid,
                          concrete_jac = true), save_everystep = false);
iter3 = iter[];
iter[] = 0;
sol4 = solve(prob_ode_brusselator_2d_sparse,
             Rosenbrock23(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid2,
                          concrete_jac = true), save_everystep = false);
iter4 = iter[];
iter[] = 0;

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test sol3.retcode === ReturnCode.Success
@test sol4.retcode === ReturnCode.Success

@test iter2 < iter1
@test iter3 < iter1
@test iter4 < iter1

sol1 = solve(prob_ode_brusselator_2d, Rodas4(linsolve = KrylovJL_GMRES()),
             save_everystep = false);
iter1 = iter[];
iter[] = 0;
sol2 = solve(prob_ode_brusselator_2d_sparse,
             Rodas4(linsolve = KrylovJL_GMRES(), precs = incompletelu, concrete_jac = true),
             save_everystep = false);
iter2 = iter[];
iter[] = 0;
sol3 = solve(prob_ode_brusselator_2d_sparse,
             Rodas4(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid,
                    concrete_jac = true), save_everystep = false);
iter3 = iter[];
iter[] = 0;
sol4 = solve(prob_ode_brusselator_2d_sparse,
             Rodas4(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid2,
                    concrete_jac = true), save_everystep = false);
iter4 = iter[];
iter[] = 0;

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test sol3.retcode === ReturnCode.Success
@test sol4.retcode === ReturnCode.Success

@test iter2 < iter1
@test iter3 < iter1
@test iter4 < iter1

sol1 = solve(prob_ode_brusselator_2d, Rodas5(linsolve = KrylovJL_GMRES()),
             save_everystep = false);
iter1 = iter[];
iter[] = 0;
sol2 = solve(prob_ode_brusselator_2d_sparse,
             Rodas5(linsolve = KrylovJL_GMRES(), precs = incompletelu, concrete_jac = true),
             save_everystep = false);
iter2 = iter[];
iter[] = 0;
sol3 = solve(prob_ode_brusselator_2d_sparse,
             Rodas5(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid,
                    concrete_jac = true), save_everystep = false);
iter3 = iter[];
iter[] = 0;
sol4 = solve(prob_ode_brusselator_2d_sparse,
             Rodas5(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid2,
                    concrete_jac = true), save_everystep = false);
iter4 = iter[];
iter[] = 0;

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test sol3.retcode === ReturnCode.Success
@test sol4.retcode === ReturnCode.Success

@test iter2 < iter1
@test iter3 < iter1
@test iter4 < iter1

sol1 = solve(prob_ode_brusselator_2d, TRBDF2(linsolve = KrylovJL_GMRES()),
             save_everystep = false);
iter1 = iter[];
iter[] = 0;
sol2 = solve(prob_ode_brusselator_2d_sparse,
             TRBDF2(linsolve = KrylovJL_GMRES(), precs = incompletelu, concrete_jac = true),
             save_everystep = false);
iter2 = iter[];
iter[] = 0;
sol3 = solve(prob_ode_brusselator_2d_sparse,
             TRBDF2(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid,
                    concrete_jac = true), save_everystep = false);
iter3 = iter[];
iter[] = 0;
sol4 = solve(prob_ode_brusselator_2d_sparse,
             TRBDF2(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid2,
                    concrete_jac = true), save_everystep = false);
iter4 = iter[];
iter[] = 0;

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test sol3.retcode === ReturnCode.Success
@test sol4.retcode === ReturnCode.Success

@test iter2 < iter1
@test iter3 < iter1
@test iter4 < iter1

sol1 = solve(prob_ode_brusselator_2d, TRBDF2(linsolve = KrylovJL_GMRES()),
             save_everystep = false);
iter1 = iter[];
iter[] = 0;
sol2 = solve(prob_ode_brusselator_2d_sparse,
             TRBDF2(linsolve = KrylovJL_GMRES(), precs = incompletelu,
                    concrete_jac = true), save_everystep = false);
iter2 = iter[];
iter[] = 0;
sol3 = solve(prob_ode_brusselator_2d_sparse,
             TRBDF2(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid,
                    concrete_jac = true), save_everystep = false);
iter3 = iter[];
iter[] = 0;
sol4 = solve(prob_ode_brusselator_2d_sparse,
             TRBDF2(linsolve = KrylovJL_GMRES(), precs = algebraicmultigrid2,
                    concrete_jac = true), save_everystep = false);
iter4 = iter[];
iter[] = 0;

@test sol1.retcode === ReturnCode.Success
@test sol2.retcode === ReturnCode.Success
@test sol3.retcode === ReturnCode.Success
@test sol4.retcode === ReturnCode.Success

@test iter2 < iter1
@test iter3 < iter1
@test iter4 < iter1
