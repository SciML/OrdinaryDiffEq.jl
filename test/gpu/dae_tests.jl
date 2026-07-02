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
using SciMLBase: ReturnCode, FullSpecialize

#=
Test goal: exercise GPU code paths for stiff/DAE solvers with mass matrices.

This is NOT a solver-quality test. We don't care whether Rosenbrock32 is a good
DAE solver; we care whether the GPU dispatch works and produces the same
result as the CPU dispatch for the same solver.

Pass criterion: GPU solution matches CPU solution (same solver) within tolerance.
CPU-vs-reference is recorded for the table but does not gate pass/fail.

The test matrix is (solver × jac_prototype × mass_matrix × linsolve). Each axis
is a simple vector at the top of the file — comment lines in/out to narrow scope.

Note on FullSpecialize: ODEFunctions are built with `FullSpecialize` to avoid a
FunctionWrappers bug where the wrapper's compiled signature is too narrow for
some solver caches, producing `llvmcall requires the compiler` errors. With
FullSpecialize the function is specialized on concrete types directly and
FunctionWrappers is bypassed.

    du[1] = -u[1]
    du[2] = -0.5*u[2]
        0 =  u[1] + u[2] - u[3]
        0 = -u[1] + u[2] - u[4]
=#

# ── Problem definition ───────────────────────────────────────────────────────

function dae!(du, u, p, t)
    return mul!(du, p, u)
end

P = [
    -1 0 0 0
    1 -0.5 0 0
    1 1 -1 0
    -1 1 0 -1
]

MASS_MATRIX = Diagonal([1, 1, 0, 0])
JAC_PROTOTYPE = sparse(map(x -> iszero(x) ? 0.0 : 1.0, P))
U0 = [1.0, 1.0, 0.5, 0.5]
TSPAN = (0.0, 5.0)

INITALG = BrownFullBasicInit()

# Shared solver kwargs. `maxiters` caps runaway integrations so a stuck
# solver fails cleanly rather than hangs indefinitely.
SOLVE_KWARGS = (; maxiters = 10_000)

# Helper: build an ODEFunction with FullSpecialize (CPU or GPU).
make_odef(; mass_matrix, jac_prototype) =
    ODEFunction(dae!; mass_matrix = mass_matrix, jac_prototype = jac_prototype)

# ── CPU reference ────────────────────────────────────────────────────────────

ODEF_CPU = make_odef(; mass_matrix = MASS_MATRIX, jac_prototype = JAC_PROTOTYPE)
PROB_CPU = ODEProblem(ODEF_CPU, U0, TSPAN, P; initializealg = INITALG)
SOL_REF = solve(PROB_CPU, Rodas5P(); SOLVE_KWARGS...)
SOL_REF_KRYLOV = solve(PROB_CPU, Rodas5P(linsolve = KrylovJL_GMRES()); SOLVE_KWARGS...)

# ── GPU problem data ─────────────────────────────────────────────────────────

U0_D = adapt(CuArray{Float64}, U0)
P_D = adapt(CuArray{Float64}, P)
MASS_MATRIX_D_DIAG = cu(MASS_MATRIX)

# ── Test matrix axes ─────────────────────────────────────────────────────────
# Comment lines to narrow scope during debugging.

# Each entry: (name, jac_prototype_for_gpu, needs_krylov_on_cpu)
JAC_VARIANTS = [
    ("none", nothing, false),
    ("CSC",  CUDA.CUSPARSE.CuSparseMatrixCSC(JAC_PROTOTYPE), true),
    ("CSR",  CUDA.CUSPARSE.CuSparseMatrixCSR(JAC_PROTOTYPE), false),
]

# Each entry: (name, gpu_mass_matrix)
MASS_VARIANTS = [
    ("diag_cu", MASS_MATRIX_D_DIAG),
    # ("dense_cu", CuArray(Matrix(MASS_MATRIX))),                           # not supported yet
    # ("sparse_cu", CUDA.CUSPARSE.CuSparseMatrixCSC(sparse(MASS_MATRIX))),  # not supported yet
]

# ── Solvers ──────────────────────────────────────────────────────────────────
# Classification reflects suitability for mass-matrix DAEs, purely informational
# for the report. Does NOT affect pass/fail — we only gate on GPU-vs-CPU match.
#
#   :suitable   — designed or known-good for index-1 mass-matrix DAEs
#   :marginal   — may lose order / show artifacts but typically runs
#   :unsuitable — included only for GPU code-path coverage

SOLVERS = [
    # ── Rosenbrock 2nd order ──
    # (Rosenbrock23,  :marginal),
    # (Rosenbrock32,  :unsuitable),   # low-accuracy, used for coverage
    (ROS2,          :marginal),
    (ROS2PR,        :suitable),
    (ROS2S,         :suitable),
    # ── Rosenbrock 3rd order ──
    # (ROS3,          :marginal),
    (ROS3PR,        :suitable),
    (ROS3PRL,       :suitable),
    (ROS3PRL2,      :suitable),
    (ROS3P,         :suitable),
    (Rodas3,        :suitable),
    # Rodas23W    — scalar indexing, needs `calculate_interpoldiff!` rework
    # Rodas3P     — scalar indexing
    (Scholz4_7,     :suitable),
    # ── Rosenbrock 4th order ──
    #=
    (ROS34PW1a,     :suitable),
    (ROS34PW1b,     :suitable),
    (ROS34PW2,      :suitable),
    (ROS34PW3,      :suitable),
    (ROS34PRw,      :suitable),
    # (RosShamp4,     :marginal),     # classical Rosenbrock, not DAE-derived
    # (Veldd4,        :marginal), # d
    (Velds4,        :marginal),
    # (GRK4T,         :marginal), # does not work with CSC/Krylov
    (GRK4A,         :marginal),
    (Ros4LStab,     :marginal),
    (Rodas4,        :suitable),
    (Rodas42,       :suitable),
    (Rodas4P,       :suitable),
    (Rodas4P2,      :suitable),
    # (ROK4a,         :suitable),
    # ── Rosenbrock 5th order ──
    (Rodas5,        :suitable),
    (Rodas5P,       :suitable),
    # (Rodas5Pe,      :suitable),
    (Rodas5Pr,      :suitable),
    # ── Rosenbrock 6th order ──
    (Rodas6P,       :suitable),
    # ── SDIRK ──
    (ImplicitEuler, :marginal),
    (Trapezoid,     :unsuitable),   # oscillates on algebraic states
    (SDIRK2,        :suitable),
    (Cash4,         :suitable),
    (Hairer4,       :suitable),
    (Hairer42,      :suitable),
    # ── BDF ──
    (ABDF2,         :suitable),
    (QNDF1,         :marginal),
    (QNDF2,         :suitable),
    # QNDF        — DeviceMemory error in LinAlg
    # (QBDF1,         :marginal),
    # (QBDF2,         :suitable),
    # QBDF        — DeviceMemory error in LinAlg
    # FBDF        — scalar indexing, needs reinitFBDF! rework
    # ── FIRK (all need `perform_step!` rework) ──
    # RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau
    #   — scalar indexing, ComplexF64 sparse unsupported
    =#
]

# ── Pass-criterion configuration ─────────────────────────────────────────────

# GPU must match CPU within this tolerance. Single knob.
GPU_MATCH_TOL = (atol = 1.0e-3, rtol = 1.0e-3)

# Known-exception list for solvers that genuinely can't meet GPU_MATCH_TOL
# (e.g. intrinsically low accuracy, non-L-stable artifacts). Entries here
# should be rare and each should be justified in a comment.
LOOSE_GPU_MATCH_TOL = Dict{Type, NamedTuple}(
    # Rosenbrock32 => (atol = Inf, rtol = Inf),  # low-accuracy method
    # Trapezoid    => (atol = 1e-2, rtol = 3e-2), # oscillation amplitude drift
)

gpu_match_tol(solver) = get(LOOSE_GPU_MATCH_TOL, solver, GPU_MATCH_TOL)

# ── Types ────────────────────────────────────────────────────────────────────

struct TestCase
    solver::Any
    solver_class::Symbol
    jac_name::String
    jac_prototype::Any
    needs_krylov_cpu::Bool
    mass_name::String
    mass_matrix::Any
end

Base.show(io::IO, c::TestCase) = print(io,
    "$(nameof(c.solver)) [jac=$(c.jac_name), mass=$(c.mass_name)]")

mutable struct TestResult
    case::TestCase
    cpu_abs::Float64      # vs reference (informational)
    cpu_rel::Float64
    gpu_abs::Float64      # vs CPU same-solver (pass/fail)
    gpu_rel::Float64
    cpu_retcode::Union{Nothing, ReturnCode.T}
    gpu_retcode::Union{Nothing, ReturnCode.T}
    cpu_error::Union{Nothing, Exception}
    gpu_error::Union{Nothing, Exception}
    status::Symbol        # :pass, :gpu_mismatch, :cpu_failed, :gpu_failed,
                          # :cpu_error, :gpu_error
end

# ── Helpers ──────────────────────────────────────────────────────────────────

function build_cases(; solvers = SOLVERS,
                       jac_variants = JAC_VARIANTS,
                       mass_variants = MASS_VARIANTS)
    cases = TestCase[]
    for (sv, cls) in solvers,
        (jn, jp, needs_krylov) in jac_variants,
        (mn, mm) in mass_variants
        push!(cases, TestCase(sv, cls, jn, jp, needs_krylov, mn, mm))
    end
    return cases
end

function maxerrs(sol_a, sol_b)
    max_abs = 0.0
    max_rel = 0.0
    for t in TSPAN[1]:0.1:TSPAN[2]
        a = Vector(sol_a(t))
        b = Vector(sol_b(t))
        diff = abs.(a .- b)
        ref = abs.(b)
        max_abs = max(max_abs, maximum(diff))
        max_rel = max(max_rel, maximum(diff ./ max.(ref, eps())))
    end
    return max_abs, max_rel
end

within_tol(abs_err, rel_err, tol) =
    abs_err < tol.atol || rel_err < tol.rtol

ok(retcode) = retcode === ReturnCode.Success || retcode === ReturnCode.Default

# ── Run a single case ────────────────────────────────────────────────────────

function run_case(case::TestCase)
    result = TestResult(case, NaN, NaN, NaN, NaN,
                        nothing, nothing, nothing, nothing, :pending)

    # CPU run ---------------------------------------------------------------
    ref = case.needs_krylov_cpu ? SOL_REF_KRYLOV : SOL_REF
    cpu_alg = case.needs_krylov_cpu ?
              case.solver(linsolve = KrylovJL_GMRES()) :
              case.solver()

    sol_cpu = nothing
    try
        sol_cpu = solve(PROB_CPU, cpu_alg; SOLVE_KWARGS...)
        result.cpu_retcode = sol_cpu.retcode
        if ok(sol_cpu.retcode)
            result.cpu_abs, result.cpu_rel = maxerrs(sol_cpu, ref)
        end
    catch e
        result.cpu_error = e
    end

    # GPU run ---------------------------------------------------------------
    sol_gpu = nothing
    try
        odef_d = make_odef(;
            mass_matrix = case.mass_matrix,
            jac_prototype = case.jac_prototype)
        prob_d = ODEProblem(odef_d, U0_D, TSPAN, P_D; initializealg = INITALG)
        sol_gpu = solve(prob_d, cpu_alg; SOLVE_KWARGS...)
        result.gpu_retcode = sol_gpu.retcode
        if ok(sol_gpu.retcode) && sol_cpu !== nothing && ok(sol_cpu.retcode)
            result.gpu_abs, result.gpu_rel = maxerrs(sol_gpu, sol_cpu)
        end
    catch e
        result.gpu_error = e
    end

    # Classify --------------------------------------------------------------
    result.status = classify(result)
    return result
end

function classify(r::TestResult)
    r.gpu_error !== nothing && return :gpu_error
    r.cpu_error !== nothing && return :cpu_error
    r.gpu_retcode !== nothing && !ok(r.gpu_retcode) && return :gpu_failed
    r.cpu_retcode !== nothing && !ok(r.cpu_retcode) && return :cpu_failed
    within_tol(r.gpu_abs, r.gpu_rel, gpu_match_tol(r.case.solver)) && return :pass
    return :gpu_mismatch
end

# ── Run & report ─────────────────────────────────────────────────────────────

function run_all(cases = build_cases(); verbose = true)
    results = TestResult[]
    Threads.@threads for (i, case) in collect(enumerate(cases))
    # for (i, case) in collect(enumerate(cases))
        verbose && @printf("[%3d/%3d] %s ... \n", i, length(cases), case)
        r = run_case(case)
        push!(results, r)
        # verbose && println(status_glyph(r.status))
    end
    return results
end

status_glyph(s) = s === :pass          ? "✓" :
                  s === :gpu_mismatch  ? "✗ gpu mismatch" :
                  s === :cpu_failed    ? "✗ cpu retcode" :
                  s === :gpu_failed    ? "✗ gpu retcode" :
                  s === :cpu_error     ? "✗ cpu error" :
                  s === :gpu_error     ? "✗ gpu error" :
                                         string(s)

function show_results(results; threshold = 1.0e-3)
    function _fmt(val)
        isnan(val) && return "   ---   "
        s = @sprintf("%.2e", val)
        return val > threshold ? "\e[31m$s\e[0m" : s
    end

    _rc(rc) = rc === nothing ? "---" : string(rc)

    label_w = 42
    println(rpad("Solver / jac / mass", label_w),
            "class       cpu_abs   cpu_rel   gpu_abs   gpu_rel   cpu_rc / gpu_rc         status")
    println("-"^140)

    for r in results
        c = r.case
        label = rpad("$(nameof(c.solver)) / $(c.jac_name) / $(c.mass_name)", label_w)
        cls = rpad(string(c.solver_class), 11)
        print(label, cls)
        print(_fmt(r.cpu_abs), "  ", _fmt(r.cpu_rel), "  ")
        print(_fmt(r.gpu_abs), "  ", _fmt(r.gpu_rel), "  ")
        print(rpad("$(_rc(r.cpu_retcode)) / $(_rc(r.gpu_retcode))", 24))
        color = r.status === :pass ? :green : :red
        printstyled(status_glyph(r.status), "\n"; color)
    end
    println("-"^140)

    # Summary
    counts = Dict{Symbol,Int}()
    for r in results
        counts[r.status] = get(counts, r.status, 0) + 1
    end
    println("\nSummary:")
    for s in (:pass, :gpu_mismatch, :cpu_failed, :gpu_failed, :cpu_error, :gpu_error)
        c = get(counts, s, 0)
        c > 0 && println("  $(rpad(string(s), 15)) $c")
    end
    return nothing
end

function show_errors(results)
    for r in results
        r.gpu_error === nothing && r.cpu_error === nothing && continue
        println("── $(r.case) ──")
        r.cpu_error !== nothing && (println("CPU error:"); showerror(stdout, r.cpu_error); println())
        r.gpu_error !== nothing && (println("GPU error:"); showerror(stdout, r.gpu_error); println())
        println()
    end
end

# ── Convenience for interactive debugging ────────────────────────────────────

"""
    debug_case(solver; jac="none", mass="diag_cu")

Run a single case interactively. Throws instead of catching so you can inspect
the stack trace.
"""
function debug_case(solver::Type; jac = "none", mass = "diag_cu")
    jv = only(filter(x -> x[1] == jac, JAC_VARIANTS))
    mv = only(filter(x -> x[1] == mass, MASS_VARIANTS))
    cls = something(SOLVERS[findfirst(s -> s[1] == solver, SOLVERS)], (solver, :unknown))[2]
    case = TestCase(solver, cls, jv[1], jv[2], jv[3], mv[1], mv[2])

    ref = case.needs_krylov_cpu ? SOL_REF_KRYLOV : SOL_REF
    cpu_alg = case.needs_krylov_cpu ?
              case.solver(linsolve = KrylovJL_GMRES()) :
              case.solver()

    @info "CPU solve"
    sol_cpu = solve(PROB_CPU, cpu_alg; SOLVE_KWARGS...)
    @info "CPU retcode" sol_cpu.retcode
    @info "GPU solve"
    odef_d = make_odef(;
        mass_matrix = case.mass_matrix,
        jac_prototype = case.jac_prototype)
    prob_d = ODEProblem(odef_d, U0_D, TSPAN, P_D; initializealg = INITALG)
    sol_gpu = solve(prob_d, case.solver(); SOLVE_KWARGS...)
    @info "GPU retcode" sol_gpu.retcode

    cpu_abs, cpu_rel = maxerrs(sol_cpu, ref)
    gpu_abs, gpu_rel = maxerrs(sol_gpu, sol_cpu)
    @info "results" cpu_abs cpu_rel gpu_abs gpu_rel
    return (; sol_cpu, sol_gpu, cpu_abs, cpu_rel, gpu_abs, gpu_rel)
end

# ── Entry point ──────────────────────────────────────────────────────────────

@testset "GPU DAE solver compatibility" begin
    global results = run_all()
    show_results(results)
    @testset "$(r.case)" for r in results
        @test r.status === :pass
    end
end

# res = debug_case(Rodas5P; jac = "CSR");
# res = debug_case(Rosenbrock23; jac = "CSC");)
