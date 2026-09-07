# FBDF Time Filter Accuracy & Work-Precision Tests
# Tests Simple FBDF_{k+1} implementation (DeCaria et al., arXiv:1810.06670v1)
#   BDF_k (order k) -> filter -> FBDF_{k+1} (order k+1)
#   y^{k+1} = y^k - η^{k+1} δ^{k+1} y^k
#
# Single flag time_filter now merges both:
#   - k==3: embedded 2-3-4 family (y2 via BDF3-Stab μ=9/125 G-stable order2,
#           y3 BDF3, y4 FBDF4 order4) from one BDF3 solve, pick max h among 2,3,4
#   - k!=3, k<=4: y_{k+1}=y_k - η δ^{k+1}, pick max h among k,k+1
# Gives larger steps: accuracy-limited -> h_{k+1}>h_k, stability-limited k==3 -> h2 may win.
#
# Measures:
#   - steps, nf, error, err/tol, avg dt (larger dt = larger step size)
#   - work-precision: error vs steps, error vs nf
#   - convergence order with fixed dt
#   - comparison: FBDF, FBDF(time_filter) [single flag includes stab+time], MOOSE234
#
# Usage:
#   julia --project=lib/OrdinaryDiffEqBDF test/fbdf_filter_accuracy.jl

using OrdinaryDiffEqBDF
using SciMLBase
using LinearAlgebra
using Printf

# ---------------------------------------------------------------------------
# Problem definitions
# ---------------------------------------------------------------------------

function exp_decay!(du,u,p,t)
    du[1] = -u[1]
end
const PROB_EXP = ODEProblem(exp_decay!, [1.0], (0.0, 10.0))
const EXACT_EXP_10 = exp(-10.0)
const EXACT_EXP_1 = exp(-1.0)

function vdp!(du,u,p,t)
    μ = p
    du[1] = u[2]
    du[2] = μ*(1 - u[1]^2)*u[2] - u[1]
end

function rober!(du,u,p,t)
    du[1] = -0.04*u[1] + 1e4*u[2]*u[3]
    du[2] =  0.04*u[1] - 1e4*u[2]*u[3] - 3e7*u[2]^2
    du[3] =  3e7*u[2]^2
end

function hires!(du,u,p,t)
    du[1] = -1.71*u[1] + 0.43*u[2] + 8.32*u[3] + 0.0007
    du[2] = 1.71*u[1] - 8.75*u[2]
    du[3] = -10.03*u[3] + 0.43*u[4] + 0.035*u[5]
    du[4] = 8.32*u[2] + 1.71*u[3] - 1.12*u[4]
    du[5] = -1.745*u[5] + 0.43*u[6] + 0.43*u[7]
    du[6] = -280.0*u[6]*u[8] + 0.69*u[4] + 1.71*u[5] - 0.43*u[6] + 0.69*u[7]
    du[7] = 280.0*u[6]*u[8] - 1.81*u[7]
    du[8] = -280.0*u[6]*u[8] + 1.81*u[7]
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

struct Metrics
    alg_name::String
    tol::Float64
    steps::Int
    nf::Int
    err::Float64
    err_over_tol::Float64
    avg_dt::Float64
    retcode::Symbol
end

function compute_metrics(sol, exact_or_ref, tol, alg_name, tspan_len)
    err = if exact_or_ref isa Number
        abs(sol.u[end][1] - exact_or_ref)
    else
        norm(sol.u[end] - exact_or_ref)
    end
    steps = length(sol.t)
    nf = sol.stats.nf
    avg_dt = tspan_len / steps
    return Metrics(alg_name, tol, steps, nf, err, err/tol, avg_dt, Symbol(sol.retcode))
end

function print_header(title)
    println("\n" * "="^80)
    println(title)
    println("="^80)
    @printf("%-12s %-10s %6s %6s %12s %10s %10s %8s\n",
        "Algorithm","tol","steps","nf","err","err/tol","avg_dt","retcode")
    println("-"^80)
end

function print_metrics(m::Metrics)
    @printf("%-12s %10.0e %6d %6d %12.3e %10.2f %10.3e %8s\n",
        m.alg_name, m.tol, m.steps, m.nf, m.err, m.err_over_tol, m.avg_dt, m.retcode)
end

# ---------------------------------------------------------------------------
# Test 1: Exp decay - exact solution, accuracy-limited, filter should give larger steps
# ---------------------------------------------------------------------------

function test_exp_decay()
    print_header("Test 1: Exp decay u'=-u, t=[0,10], exact=exp(-10), accuracy-limited")
    println("# Single time_filter flag now includes both stab (k==3->2) and time (k->k+1)")
    println("# Filter should give fewer steps (larger avg_dt) with similar err/tol<1")
    println("# At tight tol, filter should also give fewer nf (higher order offsets extra f eval)")

    for tol in [1e-4, 1e-6, 1e-8, 1e-10]
        tspan_len = 10.0
        sol1 = solve(PROB_EXP, FBDF(), reltol=tol, abstol=tol)
        m1 = compute_metrics(sol1, EXACT_EXP_10, tol, "FBDF", tspan_len)
        print_metrics(m1)

        sol2 = solve(PROB_EXP, FBDF(time_filter=true), reltol=tol, abstol=tol)
        m2 = compute_metrics(sol2, EXACT_EXP_10, tol, "FBDF+TF", tspan_len)
        print_metrics(m2)

        try
            sol3 = solve(PROB_EXP, MOOSE234(), reltol=tol, abstol=tol)
            m3 = compute_metrics(sol3, EXACT_EXP_10, tol, "MOOSE234", tspan_len)
            print_metrics(m3)
        catch e
            @printf("%-12s %10.0e  -- MOOSE234 not available or failed: %s\n", "MOOSE234", tol, e)
        end

        # Larger step check: filter should have larger avg_dt (fewer steps) for smooth problem
        if m2.steps < m1.steps
            println("  ✓ TF larger steps: $(m1.steps)→$(m2.steps) steps, avg_dt $(m1.avg_dt)→$(m2.avg_dt) ($(m2.avg_dt/m1.avg_dt)x)")
        else
            println("  ✗ TF did not reduce steps: $(m1.steps) vs $(m2.steps)")
        end
        println()
    end
end

# ---------------------------------------------------------------------------
# Test 2: Van der Pol μ=10 - moderately stiff, mixed accuracy/stability
# ---------------------------------------------------------------------------

function test_vdp()
    print_header("Test 2: Van der Pol μ=10, t=[0,20], ref=FBDF 1e-12")
    println("# Single time_filter includes k==3 -> 2 (G-stable) and k->k+1, picks max h among 2,3,4 when k==3")
    u0=[2.0,0.0]; tspan=(0.0,20.0); p=10.0
    prob=ODEProblem(vdp!,u0,tspan,p)
    tspan_len = 20.0
    ref = solve(prob, FBDF(), reltol=1e-12, abstol=1e-12)
    ref_end = ref.u[end]

    for tol in [1e-4,1e-6,1e-8]
        sol1=solve(prob, FBDF(), reltol=tol, abstol=tol)
        m1=compute_metrics(sol1, ref_end, tol, "FBDF", tspan_len)
        print_metrics(m1)

        sol2=solve(prob, FBDF(time_filter=true), reltol=tol, abstol=tol)
        m2=compute_metrics(sol2, ref_end, tol, "FBDF+TF", tspan_len)
        print_metrics(m2)

        println("  TF steps ratio: $(m2.steps/m1.steps), nf ratio: $(m2.nf/m1.nf), err ratio: $(m2.err/m1.err)")
        println()
    end
end

# ---------------------------------------------------------------------------
# Test 3: ROBER - very stiff, stability-limited, filter may hurt
# ---------------------------------------------------------------------------

function test_rober()
    print_header("Test 3: ROBER t=[0,1e4], ref=FBDF 1e-12, stability-limited")
    println("# Expect filter to increase steps/nf because FBDF4 has smaller stability region than BDF3")
    u0=[1.0,0.0,0.0]; tspan=(0.0,1e4)
    prob=ODEProblem(rober!,u0,tspan)
    tspan_len=1e4
    ref = solve(prob, FBDF(), reltol=1e-12, abstol=1e-12)
    ref_end = ref.u[end]

    for tol in [1e-6,1e-8,1e-10]
        sol1=solve(prob, FBDF(), reltol=tol, abstol=tol)
        m1=compute_metrics(sol1, ref_end, tol, "FBDF", tspan_len)
        print_metrics(m1)

        sol2=solve(prob, FBDF(time_filter=true), reltol=tol, abstol=tol)
        m2=compute_metrics(sol2, ref_end, tol, "FBDF+TF", tspan_len)
        print_metrics(m2)

        println("  TF/Plain steps=$(m2.steps/m1.steps), nf=$(m2.nf/m1.nf)")
        println()
    end
end

# ---------------------------------------------------------------------------
# Test 4: HIRES - 8 eq, mildly stiff, filter may help at tight tol
# ---------------------------------------------------------------------------

function test_hires()
    print_header("Test 4: HIRES t=[0,321.8122], ref=FBDF 1e-11")
    u0=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0007]
    tspan=(0.0,321.8122)
    prob=ODEProblem(hires!,u0,tspan)
    tspan_len=321.8122
    ref = solve(prob, FBDF(), reltol=1e-11, abstol=1e-11)
    ref_end = ref.u[end]

    for tol in [1e-6,1e-8]
        sol1=solve(prob, FBDF(), reltol=tol, abstol=tol)
        m1=compute_metrics(sol1, ref_end, tol, "FBDF", tspan_len)
        print_metrics(m1)

        sol2=solve(prob, FBDF(time_filter=true), reltol=tol, abstol=tol)
        m2=compute_metrics(sol2, ref_end, tol, "FBDF+TF", tspan_len)
        print_metrics(m2)

        println("  TF should give fewer steps at tight tol if accuracy-limited")
        println()
    end
end

# ---------------------------------------------------------------------------
# Test 5: Fixed dt convergence order - verifies filter raises order
# ---------------------------------------------------------------------------

function test_convergence_order()
    print_header("Test 5: Fixed dt convergence order, u'=-u, t=[0,1], exact=exp(-1)")
    println("# With adaptive=false, dt fixed, order should be ~k for FBDF, ~k+1 for FBDF+TF when k<max_order")
    println("# For k=1 (BDF1) filter should give order2, k=2->3, etc. Ratio err(dt)/err(dt/2) ~ 2^{p+1}")

    f(u,p,t) = -u
    exact = exp(-1)
    tspan=(0.0,1.0)
    u0=[1.0]

    for max_order in [2,3,4]
        println("\n-- max_order=$max_order --")
        prev_err1=nothing; prev_err2=nothing
        for dt in [0.2,0.1,0.05,0.025]
            prob=ODEProblem(f,u0,tspan)
            sol1=solve(prob, FBDF(max_order=Val(max_order)), dt=dt, adaptive=false)
            err1=abs(sol1.u[end][1]-exact)
            sol2=solve(prob, FBDF(max_order=Val(max_order),time_filter=true), dt=dt, adaptive=false)
            err2=abs(sol2.u[end][1]-exact)
            if prev_err1 !== nothing
                order1 = log(prev_err1/err1)/log(2)
                order2 = log(prev_err2/err2)/log(2)
                @printf("dt=%5.3f FBDF err=%10.3e order~%4.2f | TF err=%10.3e order~%4.2f\n", dt, err1, order1, err2, order2)
            else
                @printf("dt=%5.3f FBDF err=%10.3e | TF err=%10.3e\n", dt, err1, err2)
            end
            prev_err1=err1; prev_err2=err2
        end
    end
end

# ---------------------------------------------------------------------------
# Test 6: Larger step size verification - avg_dt should be larger with filter for smooth
# ---------------------------------------------------------------------------

function test_larger_step()
    print_header("Test 6: Larger step size verification")
    println("# For smooth accuracy-limited problems, time filter should give larger avg_dt = tspan/steps")
    prob = ODEProblem((u,p,t)-> -u, [1.0], (0.0, 10.0))
    tspan_len=10.0
    for tol in [1e-6,1e-8,1e-10]
        sol1=solve(prob, FBDF(), reltol=tol, abstol=tol)
        sol2=solve(prob, FBDF(time_filter=true), reltol=tol, abstol=tol)
        avg1=tspan_len/length(sol1.t)
        avg2=tspan_len/length(sol2.t)
        println("tol=$tol FBDF avg_dt=$avg1 steps=$(length(sol1.t)) | TF avg_dt=$avg2 steps=$(length(sol2.t)) ratio=$(avg2/avg1)")
        if avg2 > avg1
            println("  ✓ TF larger avg_dt as expected (filter raises order, smaller LTE)")
        else
            println("  ✗ TF not larger - may be stability-limited or startup overhead")
        end
    end
end

# ---------------------------------------------------------------------------
# Test 7: Mass matrix / DAE - filter should not break DAE
# ---------------------------------------------------------------------------

function test_mass_matrix()
    print_header("Test 7: Mass matrix DAE, filter should handle or fallback gracefully")
    # Simple DAE: M*u' = f(u), M = [1 0; 0 0], algebraic constraint
    # du1/dt = -u1 + u2, 0 = u1 - u2  => u1=u2, du1/dt=0 => u1=const
    # Use Rosenbrock-like? For FBDF with mass matrix, filter should be skipped or use heuristic
    # Here we just test that it doesn't error
    function f!(du,u,p,t)
        du[1] = -u[1] + u[2]
        du[2] = u[1] - u[2]
    end
    M = [1.0 0.0; 0.0 0.0]
    prob = ODEProblem(ODEFunction(f!, mass_matrix=M), [1.0,1.0], (0.0,1.0))
    try
        sol1=solve(prob, FBDF(), reltol=1e-6, abstol=1e-6)
        println("FBDF mass_matrix: retcode=$(sol1.retcode) steps=$(length(sol1.t))")
        sol2=solve(prob, FBDF(time_filter=true), reltol=1e-6, abstol=1e-6)
        println("FBDF+TF mass_matrix: retcode=$(sol2.retcode) steps=$(length(sol2.t))")
        println("  ✓ No crash with mass matrix")
    catch e
        println("  ✗ Failed with mass matrix: $e")
    end
end

# ---------------------------------------------------------------------------
# Run all
# ---------------------------------------------------------------------------

function run_all()
    test_exp_decay()
    test_vdp()
    test_rober()
    test_hires()
    test_convergence_order()
    test_larger_step()
    test_mass_matrix()
    println("\n" * "="^80)
    println("All tests done. Check err/tol<1 for accuracy, steps ratio for larger step size.")
    println("Filter benefits when accuracy-limited (smooth, tight tol), not when stability-limited.")
    println("="^80)
end

run_all()
