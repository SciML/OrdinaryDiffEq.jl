using OrdinaryDiffEqAdamsBashforthMoulton, DiffEqBase, Test

# Callback `resize!(integrator, n)` must grow/shrink every Adams cache buffer that
# participates in the step (rate history, nested BS3/RK4 startup caches). The
# `@cache`-generated `full_cache` only walks `uType`/`rateType`/`uNoUnitsType`
# fields, so VCAB* ϕ-history (`coefType`) and nested caches need
# `resize_non_user_cache!`. Fixed-step AB*/ABM* already expose their rate buffers
# via `full_cache`; they are included here as a regression fence.

function f_resize!(du, u, p, t)
    @inbounds for i in eachindex(u)
        du[i] = -u[i]
    end
    return nothing
end

const u0_resize = [1.0, 2.0, 3.0]
const t_resize = 0.5
const tspan_resize = (0.0, 1.0)

exact_grow(t) = t < t_resize ? [exp(-t), 2exp(-t), 3exp(-t)] :
    [exp(-t), 2exp(-t), 3exp(-t), exp(-(t - t_resize))]
exact_shrink(t) = t < t_resize ? [exp(-t), 2exp(-t), 3exp(-t)] :
    [exp(-t), 2exp(-t)]

function solve_with_resize(alg, n_new; kwargs...)
    function affect!(integrator)
        resize!(integrator, n_new)
        n_new > length(u0_resize) && (integrator.u[n_new] = 1.0)
        return nothing
    end
    cb = DiscreteCallback((u, t, integrator) -> t == t_resize, affect!)
    return solve(
        ODEProblem(f_resize!, copy(u0_resize), tspan_resize), alg;
        callback = cb, tstops = [t_resize], kwargs...
    )
end

function max_post_resize_error(sol, exact)
    err = 0.0
    for (t, u) in zip(sol.t, sol.u)
        t == t_resize && continue
        e = exact(t)
        length(u) == length(e) || continue
        err = max(err, maximum(abs, u .- e))
    end
    return err
end

const FIXED_ABM_TYPES = Union{AB3, AB4, AB5, ABM32, ABM43, ABM54}

@testset "resize! grow/shrink: $(nameof(typeof(alg)))" for alg in (
        AB3(), AB4(), AB5(), ABM32(), ABM43(), ABM54(),
        VCAB3(), VCAB4(), VCAB5(), VCABM3(), VCABM4(), VCABM5(), VCABM(),
    )
    fixed = alg isa FIXED_ABM_TYPES
    kwargs = fixed ? (adaptive = false, dt = 0.01) :
        (abstol = 1.0e-10, reltol = 1.0e-10)
    # Fixed-step AB3/ABM32 startup is O(dt^3); adaptive VCAB* meet tighter tol.
    atol = fixed ? 1.0e-5 : 1.0e-6

    @testset "grow" begin
        sol = solve_with_resize(alg, 4; kwargs...)
        @test SciMLBase.successful_retcode(sol)
        @test length(sol.u[end]) == 4
        @test max_post_resize_error(sol, exact_grow) < atol
    end

    @testset "shrink" begin
        sol = solve_with_resize(alg, 2; kwargs...)
        @test SciMLBase.successful_retcode(sol)
        @test length(sol.u[end]) == 2
        @test max_post_resize_error(sol, exact_shrink) < atol
    end
end
