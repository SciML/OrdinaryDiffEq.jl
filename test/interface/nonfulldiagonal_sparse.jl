using OrdinaryDiffEq, DiffEqOperators, SparseArrays, LinearSolve
using UnPack
using ComponentArrays
using Symbolics

function enclosethetimedifferential(parameters::NamedTuple)::Function
    @info "Enclosing the time differential"

    @unpack Δr, r_space, countorderapprox = parameters.compute
    countdiscretizationsteps = length(r_space)

    ∇ = CenteredDifference(1, countorderapprox, Δr, countdiscretizationsteps)
    Δ = CenteredDifference(2, countorderapprox, Δr, countdiscretizationsteps)
    # To concretize as arrays
    tmpbc = NeumannBC((3.14, 1.0), Δr)
    Δ, _ = Array(Δ * tmpbc)
    ∇, _ = Array(∇ * tmpbc)
    bc_x = zeros(Real, countdiscretizationsteps)
    bc_xx = zeros(Real, countdiscretizationsteps)

    function timedifferentialclosure!(du, u, p, t)
        @unpack (α, D, v, k_p, V_c, Q_l, Q_r, V_b,
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
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> dudt(du, u, parameters.prior, 0.0),
                                           du0,
                                           odeprob.u0);
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
