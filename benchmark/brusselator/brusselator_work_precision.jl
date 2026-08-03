

using OrdinaryDiffEq, OrdinaryDiffEqBDF, OrdinaryDiffEqExponentialRK, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK
using LinearAlgebra, SparseArrays, SciMLBase, Printf

# MatrixOperator import - location moved in latest master, try all
try
    using DiffEqBase: MatrixOperator
catch
    using SciMLOperators: MatrixOperator
end


const NS = [16, 32]  # Grid sizes: 16=>512, 32=>2048, 64=>8192 ODEs (add 64 after reviewing 32)
const TOLS = [1e-4, 1e-5, 1e-6]  # Adaptive implicit tol sweep
const KRYLOV_M = 30

const EXP_ALGS = [
    ("ETDRK2_m30", () -> ETDRK2(krylov=true, m=KRYLOV_M)),
    ("ETDRK4_m30", () -> ETDRK4(krylov=true, m=KRYLOV_M)),
    ("HochOst4_m30", () -> HochOst4(krylov=true, m=KRYLOV_M)),
    ("LawsonEuler_m30", () -> LawsonEuler(krylov=true, m=KRYLOV_M)),
]

const IMPLICIT_ALGS = [
    ("FBDF", () -> FBDF()),
    ("Rodas5P", () -> Rodas5P()),
    ("KenCarp4", () -> KenCarp4()),
    ("TRBDF2", () -> TRBDF2()),
]

function estimate_builds(name::String)
    if occursin("ETDRK2", name); 2
    elseif occursin("ETDRK3", name); 4
    elseif occursin("ETDRK4", name); 5
    elseif occursin("HochOst4", name); 9
    else; 1
    end
end

# =============================================================================
# Problem definition - matrix operator version (fixes size(Function) bug)
# =============================================================================
brusselator_f(x, y, t) = (((x - 0.3)^2 + (y - 0.6)^2) <= 0.1^2) * (t >= 1.1) * 5.0

function init_brusselator_2d(N::Int)
    xyd = range(0, stop=1, length=N)
    u = zeros(N, N, 2)
    for I in CartesianIndices((N, N))
        x = xyd[I[1]]; y = xyd[I[2]]
        u[I, 1] = 22 * (y * (1 - y))^(3 / 2)
        u[I, 2] = 27 * (x * (1 - x))^(3 / 2)
    end
    return u
end

function build_brusselator_laplacian(N::Int, alpha_scaled::Float64)
    n = N * N
    L = spzeros(n, n)
    for j in 1:N, i in 1:N
        idx = (j-1)*N + i
        ip1 = i == N ? 1 : i+1
        im1 = i == 1 ? N : i-1
        jp1 = j == N ? 1 : j+1
        jm1 = j == 1 ? N : j-1
        idx_ip1 = (j-1)*N + ip1
        idx_im1 = (j-1)*N + im1
        idx_jp1 = (jp1-1)*N + i
        idx_jm1 = (jm1-1)*N + i
        L[idx, idx] = -4 * alpha_scaled
        L[idx, idx_ip1] = alpha_scaled
        L[idx, idx_im1] = alpha_scaled
        L[idx, idx_jp1] = alpha_scaled
        L[idx, idx_jm1] = alpha_scaled
    end
    C = spzeros(2*n, 2*n)
    C[1:n, 1:n] = L
    C[n+1:end, n+1:end] = L
    return C
end

function brusselator_2d_nonlin_vec!(du, u, p, t)
    N = p[5]
    A, B = p[1], p[2]
    u3 = reshape(u, N, N, 2)
    du3 = reshape(du, N, N, 2)
    @inbounds for I in CartesianIndices((N, N))
        i, j = Tuple(I)
        x = (i-1)/N; y = (j-1)/N
        du3[i, j, 1] = B + u3[i, j, 1]^2 * u3[i, j, 2] - (A + 1) * u3[i, j, 1] + brusselator_f(x, y, t)
        du3[i, j, 2] = A * u3[i, j, 1] - u3[i, j, 1]^2 * u3[i, j, 2]
    end
end

make_p(N::Int; A=3.4, B=1.0, α=10.0) = (A, B, α, 1.0/N, N)

function make_split_prob(N::Int; tspan=(0.0, 11.5))
    p = make_p(N)
    αs = p[3] / p[4]^2
    L = build_brusselator_laplacian(N, αs)
    u0_vec = vec(init_brusselator_2d(N))
    f1 = MatrixOperator(L)
    f2 = (du, u, par, t) -> brusselator_2d_nonlin_vec!(du, u, par, t)
    return SplitODEProblem(f1, f2, u0_vec, tspan, p)
end

function make_full_prob(N::Int; tspan=(0.0, 11.5))
    p = make_p(N)
    u0_vec = vec(init_brusselator_2d(N))
    function full!(du, u, par, t)
        Nloc = par[5]
        α, dx = par[3], par[4]
        αs = α / dx^2
        u3 = reshape(u, Nloc, Nloc, 2)
        du3 = reshape(du, Nloc, Nloc, 2)
        @inbounds for k in 1:2, j in 1:Nloc, i in 1:Nloc
            ip1 = i == Nloc ? 1 : i+1; im1 = i == 1 ? Nloc : i-1
            jp1 = j == Nloc ? 1 : j+1; jm1 = j == 1 ? Nloc : j-1
            du3[i, j, k] = αs * (u3[im1, j, k] + u3[ip1, j, k] + u3[i, jp1, j, k] + u3[i, jm1, j, k] - 4 * u3[i, j, k])
        end
        tmp = similar(du)
        brusselator_2d_nonlin_vec!(tmp, u, par, t)
        du .+= tmp
    end
    return ODEProblem(full!, u0_vec, tspan, p)
end

# --- Warmup to avoid JIT compilation inflating timing (e.g., 7.1s vs 0.08s for same work) ---
function warmup_all(prob_full, prob_split)
    @printf("\n--- Warmup: dummy solves to precompile (JIT) ---\n")
    for (name, factory) in IMPLICIT_ALGS
        try solve(prob_full, factory(), abstol=1e-4, reltol=1e-4, maxiters=10, save_everystep=false) catch; end
    end
    for (name, factory) in EXP_ALGS
        try solve(prob_split, factory(), dt=0.05, maxiters=10, save_everystep=false) catch; end
    end
    @printf("Warmup done\n")
end

function run_for_N(N::Int)
    @printf("\n==================== N=%d (%d ODEs) ====================\n", N, 2*N^2)
    prob_full = make_full_prob(N)
    prob_split = make_split_prob(N)

    # Reference at tight tolerance for error measurement
    @printf("Reference FBDF 1e-12...\n")
    sol_ref = solve(prob_full, FBDF(), abstol=1e-12, reltol=1e-12, maxiters=Int(5e6), save_everystep=false)
    u_ref = sol_ref.u[end]
    @printf("Ref steps=%d\n", sol_ref.stats.naccept)

    warmup_all(prob_full, prob_split)

    results=[]
    for tol in TOLS
        for (name, factory) in IMPLICIT_ALGS
            alg=factory()
            try
                t_s=time(); sol=solve(prob_full, alg, abstol=tol, reltol=tol, maxiters=Int(5e6), save_everystep=false); t_e=time()-t_s
                err=norm(sol.u[end]-u_ref)/max(1.0,norm(u_ref))
                @printf("  %-20s tol=%.0e err=%.2e time=%6.3f nf=%d acc=%d rej=%d\n", name, tol, err, t_e, sol.stats.nf+sol.stats.nf2, sol.stats.naccept, sol.stats.nreject)
                push!(results,(N=N,tol=tol,alg=name,err=err,time=t_e,nf=sol.stats.nf+sol.stats.nf2,naccept=sol.stats.naccept,nreject=sol.stats.nreject,njacs=sol.stats.njacs,nw=sol.stats.nw,builds=0,m=0,est_matvecs=0))
            catch e; @printf("  %-20s EX %s\n", name, e); end
        end
        for (name, factory) in EXP_ALGS
            alg=factory()
            try
                t_s=time(); sol=solve(prob_split, alg, dt=0.05, abstol=tol, reltol=tol, maxiters=Int(5e6), save_everystep=false); t_e=time()-t_s
                err=norm(sol.u[end]-u_ref)/max(1.0,norm(u_ref))
                m_val=hasproperty(alg,:m) ? alg.m : 0; builds=estimate_builds(name)
                @printf("  %-20s tol=%.0e err=%.2e time=%6.3f nf=%d acc=%d est_matvec=%d\n", name, tol, err, t_e, sol.stats.nf+sol.stats.nf2, sol.stats.naccept, builds*m_val*sol.stats.naccept)
                push!(results,(N=N,tol=tol,alg=name,err=err,time=t_e,nf=sol.stats.nf+sol.stats.nf2,naccept=sol.stats.naccept,nreject=sol.stats.nreject,njacs=sol.stats.njacs,nw=sol.stats.nw,builds=builds,m=m_val,est_matvecs=builds*m_val*sol.stats.naccept))
            catch e; @printf("  %-20s EX %s\n", name, e); end
        end
    end
    return results
end

all_results=[]
for N in NS
    append!(all_results, run_for_N(N))
end

csv_path=joinpath(@__DIR__, "brusselator_benchmark_results.csv")
open(csv_path,"w") do io
    println(io,"N,tol,alg,err,time,nf,naccept,nreject,njacs,nw,builds,m,est_matvecs")
    for r in all_results
        println(io,"$(r.N),$(r.tol),$(r.alg),$(r.err),$(r.time),$(r.nf),$(r.naccept),$(r.nreject),$(r.njacs),$(r.nw),$(r.builds),$(r.m),$(r.est_matvecs)")
    end
end
println("\nCSV saved: $csv_path")
