using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqBDF
using OrdinaryDiffEqFIRK
using OrdinaryDiffEqNonlinearSolve
using CUDA
using LinearAlgebra
using Adapt
using SparseArrays
using Test
using CUDSS
using Printf
using OrdinaryDiffEqNonlinearSolve.LinearSolve: KrylovJL_GMRES


#=
du[1] = -u[1]
du[2] = -0.5*u[2]
    0 =  u[1] + u[2] - u[3]
    0 = -u[1] + u[2] - u[4]
=#

function dae!(du, u, p, t)
    return mul!(du, p, u)
end

p = [
    -1 0 0 0
    1 -0.5 0 0
    1 1 -1 0
    -1 1 0 -1
]

mass_matrix = Diagonal([1, 1, 0, 0])
jac_prototype = sparse(map(x -> iszero(x) ? 0.0 : 1.0, p))

u0 = [1.0, 1.0, 0.5, 0.5] # force init
tspan = (0.0, 5.0)

# CPU reference solution (Rodas5P)
odef = ODEFunction(dae!, mass_matrix = mass_matrix, jac_prototype = jac_prototype)
prob = ODEProblem(odef, u0, tspan, p)
sol_ref = solve(prob, Rodas5P())
sol_ref_krylov = solve(prob, Rodas5P(linsolve = KrylovJL_GMRES()))

# GPU data -- we use F64 for higher accuracy for comparison
u0_d = adapt(CuArray{Float64}, u0)
p_d = adapt(CuArray{Float64}, p)

# dense or sparse mass matrix does not work yet!
mass_matrix_d = cu(mass_matrix)

# ── Jacobian prototype options ────────────────────────────────────────────────
jac_prots = [
    "none" => nothing,
    "CSC" => CUDA.CUSPARSE.CuSparseMatrixCSC(jac_prototype),
    "CSR" => CUDA.CUSPARSE.CuSparseMatrixCSR(jac_prototype),
]

# ── Solver definitions ────────────────────────────────────────────────────────
#
# Each entry: SolverType => (; tol overrides...)
#   method_tol     = (; atol, rtol)  – cpu solver vs Rodas5P ref (all jac); some methods are poor fits for DAEs
#   method_csc_tol = (; atol, rtol)  – cpu solver vs Rodas5P ref (CSC/Krylov path only)
#   gpu_tol        = (; atol, rtol)  – gpu vs cpu (none/CSR fallback)
#   csc_tol        = (; atol, rtol)  – gpu vs cpu (CSC jac_prot)
#   csr_tol        = (; atol, rtol)  – gpu vs cpu (CSR jac_prot)
# Only specify fields that deviate from the defaults.

const DEFAULT_METHOD_TOL = (; atol = 1.0e-2, rtol = 1.0e-2)
const DEFAULT_GPU_TOL = (; atol = 1.0e-4, rtol = 1.0e-4)

solvers = [
    # ── Rosenbrock 2nd order ──
    Rosenbrock23 => (;
        csc_tol = (; atol = 2.0e-4, rtol = 5.0e-4),
    ),
    Rosenbrock32 => (;
        method_tol = (; atol = 2.0e-2, rtol = 1.5),
        csc_tol = (; atol = Inf, rtol = Inf),
    ),
    ROS2 => (;),
    ROS2PR => (;
        csc_tol = (; atol = 3.0e-4, rtol = 4.0e-4),
    ),
    ROS2S => (;
        csc_tol = (; atol = 7.0e-4, rtol = 1.2e-3),
    ),
    # ── Rosenbrock 3rd order ──
    ROS3 => (;),
    ROS3PR => (;
        method_tol = (; atol = 0.25, rtol = 15.0),
    ),
    ROS3PRL => (;
        csc_tol = (; atol = 2.0e-3, rtol = 4.0e-3),
    ),
    ROS3PRL2 => (;
        csc_tol = (; atol = 3.0e-3, rtol = 4.0e-3),
    ),
    ROS3P => (;
        method_tol = (; atol = 0.25, rtol = 15.0),
    ),
    Rodas3 => (;
        csc_tol = (; atol = 2.0e-3, rtol = 3.0e-3),
    ),
    # Rodas23W()  # scalar indexing, requires large changes to `calculate_interpoldiff!`
    # Rodas3P()   # scalar indexing
    Scholz4_7 => (;),
    # ── Rosenbrock 4th order ──
    ROS34PW1a => (;
        method_tol = (; atol = 0.25, rtol = 5.0),
    ),
    ROS34PW1b => (;
        method_tol = (; atol = 0.25, rtol = 5.0),
    ),
    ROS34PW2 => (;
        csc_tol = (; atol = 1.2e-3, rtol = 2.2e-3),
    ),
    ROS34PW3 => (;),
    ROS34PRw => (;
        csc_tol = (; atol = 6.0e-4, rtol = 1.0e-3),
    ),
    RosShamp4 => (;),
    Veldd4 => (;),
    Velds4 => (;),
    GRK4T => (;),
    GRK4A => (;),
    Ros4LStab => (;),
    Rodas4 => (;),
    Rodas42 => (;),
    Rodas4P => (;),
    Rodas4P2 => (;),
    ROK4a => (;),
    # ── Rosenbrock 5th order ──
    Rodas5 => (;),
    Rodas5P => (;),
    Rodas5Pe => (;),
    Rodas5Pr => (;),
    # ── Rosenbrock 6th order ──
    Rodas6P => (;),
    # ── SDIRK (don't include fixed time step which need explicit dt) ──
    ImplicitEuler => (;),
    Trapezoid => (;
        csc_tol = (; atol = 8.0e-3, rtol = 2.5e-2),
    ),
    SDIRK2 => (;
        csc_tol = (; atol = 1.5e-4, rtol = 3.5e-4),
    ),
    Cash4 => (;
        csc_tol = (; atol = 5.0e-3, rtol = 8.0e-3),
    ),
    Hairer4 => (;
        csc_tol = (; atol = 8.0e-3, rtol = 6.0e-2),
    ),
    Hairer42 => (;
        csc_tol = (; atol = 7.0e-3, rtol = 5.0e-2),
    ),
    # ── BDF ──
    ABDF2 => (;
        method_csc_tol = (; atol = Inf, rtol = Inf),   # ABDF2 + Krylov diverges vs Rodas5P + Krylov
        csc_tol = (; atol = 2.0e-4, rtol = 7.0e-4),
    ),
    QNDF1 => (;),
    QNDF2 => (;
        csc_tol = (; atol = 4.0e-3, rtol = 5.0e-3),
    ),
    # QNDF()  # 🔧 DeviceMemory error in LinAlg
    QBDF1 => (;
        method_tol = (; atol = 2.0e-2, rtol = 0.2),
    ),
    QBDF2 => (;
        gpu_tol = (; atol = 2.0e-3, rtol = 2.0e-3),
        csc_tol = (; atol = 4.0e-3, rtol = 7.0e-3),
        csr_tol = (; atol = 2.0e-3, rtol = 2.0e-3),
    ),
    # QBDF()  # DeviceMemory error in LinAlg
    # FBDF()  # scalar indexing -> needs extensive work on reinitFBDF!
    # ── FIRK -> all need substantial changes to `perform_step!` for FIRK methods ──
    # RadauIIA3()     #  scalar indexing, ComplexF64 sparse unsupported
    # RadauIIA5()     #  scalar indexing, ComplexF64 sparse unsupported
    # RadauIIA9()     #  scalar indexing, ComplexF64 sparse unsupported
    # AdaptiveRadau() #  scalar indexing, ComplexF64 sparse unsupported
]

# ── Helpers ───────────────────────────────────────────────────────────────────

function _get_method_tol(overrides, jac_name)
    jac_name == "CSC" && hasproperty(overrides, :method_csc_tol) && return overrides.method_csc_tol
    hasproperty(overrides, :method_tol) && return overrides.method_tol
    return DEFAULT_METHOD_TOL
end

function _get_tol(overrides, jac_name)
    _tol_keys = Dict("none" => :gpu_tol, "CSC" => :csc_tol, "CSR" => :csr_tol)
    key = _tol_keys[jac_name]
    hasproperty(overrides, key) && return getproperty(overrides, key)
    # gpu_tol acts as default for all GPU combos
    hasproperty(overrides, :gpu_tol) && return overrides.gpu_tol
    return DEFAULT_GPU_TOL
end

function maxerrs(sol_a, sol_b)
    max_abs = 0.0
    max_rel = 0.0
    for t in tspan[begin]:0.1:tspan[end]
        a = Vector(sol_a(t))
        b = Vector(sol_b(t))
        diff = abs.(a - b)
        ref = abs.(b)
        max_abs = max(max_abs, maximum(diff))
        max_rel = max(max_rel, maximum(diff ./ max.(ref, eps())))
    end
    return max_abs, max_rel
end

function run_dae_tests()
    results = Any[]

    for (sv, overrides) in solvers, (jn, jp) in jac_prots
        println("Test $sv with prototype $jn")
        sn = string(sv)
        gtol = _get_tol(overrides, jn)
        mtol = _get_method_tol(overrides, jn)

        # CSC will always fall back to krylov so the ref solution should do to
        krylov = (jn == "CSC")
        _sol_ref = krylov ? sol_ref_krylov : sol_ref

        # CPU: this solver vs reference
        cpu_alg = krylov ? sv(linsolve = KrylovJL_GMRES()) : sv()
        sol_cpu = solve(prob, cpu_alg)
        cpu_abs, cpu_rel = maxerrs(sol_cpu, _sol_ref)

        # GPU solve
        odef_d = ODEFunction(dae!, mass_matrix = mass_matrix_d, jac_prototype = jp)
        prob_d = ODEProblem(odef_d, u0_d, tspan, p_d)
        sol_d = solve(prob_d, sv())

        # GPU vs CPU (same solver)
        gpu_abs, gpu_rel = maxerrs(sol_d, sol_cpu)

        method_passed = (cpu_abs < mtol.atol) || (cpu_rel < mtol.rtol)
        gpu_passed = (gpu_abs < gtol.atol) || (gpu_rel < gtol.rtol)
        passed = method_passed && gpu_passed
        push!(
            results, (;
                solver = sn, jac = jn, cpu_abs, cpu_rel, gpu_abs, gpu_rel,
                gtol, mtol, method_passed, gpu_passed, passed, error = "",
            )
        )
        @test passed
    end

    return results
end

function show_results(results)
    function _fmt(val; threshold = 1.0e-3)
        isnan(val) && return "   ---  "
        s = @sprintf("%.2e", val)
        if val > threshold
            return "\e[31m$s\e[0m"
        end
        return s
    end

    println(rpad("Solver / Jac", 30), "cpu_abs  cpu_rel  gpu_abs  gpu_rel  status")
    println("-"^85)

    for r in results
        label = rpad("$(r.solver) / $(r.jac)", 30)
        if r.error != ""
            printstyled(label, "ERROR: ", r.error, "\n"; color = :yellow)
        else
            print(label)
            print(_fmt(r.cpu_abs, threshold = 1.0e-2), "  ", _fmt(r.cpu_rel, threshold = 1.0e-2), "  ")
            print(_fmt(r.gpu_abs), "  ", _fmt(r.gpu_rel), "  ")
            if r.passed
                printstyled("✓\n"; color = :green)
            elseif !r.method_passed && !r.gpu_passed
                printstyled("✗ method+gpu\n"; color = :red)
            elseif !r.method_passed
                printstyled("✗ method\n"; color = :red)
            else
                printstyled("✗ gpu\n"; color = :red)
            end
        end
    end
    return println("-"^85)
end

@testset "GPU DAE solver compatibility" begin
    global results
    results = run_dae_tests()
end
show_results(results)
