using OrdinaryDiffEqAdamsBashforthMoulton, DiffEqBase, LinearAlgebra, Test

# Coupled linear resize: ϕ history / nested BS3/RK4 caches must track state size.

const A_resize = [
    -1.0 0.5 0.0 0.3;
    0.2 -2.0 0.4 0.0;
    0.0 0.3 -1.5 0.6;
    0.7 0.0 0.5 -1.0
]

function f_resize!(du, u, p, t)
    n = length(u)
    mul!(du, view(A_resize, 1:n, 1:n), u)
    return nothing
end

const u0_resize = [1.0, 2.0, 3.0]
const t_resize = 0.5
const tspan_resize = (0.0, 2.0)
const u_new_grow = 1.0

function exact_final(n_new)
    u = exp(A_resize[1:3, 1:3] * t_resize) * u0_resize
    v = n_new > 3 ? vcat(u, u_new_grow) : u[1:n_new]
    return exp(A_resize[1:n_new, 1:n_new] * (tspan_resize[2] - t_resize)) * v
end

function exact_at(t, n_new)
    u = exp(A_resize[1:3, 1:3] * min(t, t_resize)) * u0_resize
    t <= t_resize && return u
    v = n_new > 3 ? vcat(u, u_new_grow) : u[1:n_new]
    return exp(A_resize[1:n_new, 1:n_new] * (t - t_resize)) * v
end

function solve_with_resize(alg, n_new; modified = true, kwargs...)
    function affect!(integrator)
        resize!(integrator, n_new)
        n_new > length(u0_resize) && (integrator.u[n_new] = u_new_grow)
        modified || derivative_discontinuity!(integrator, false)
        return nothing
    end
    cb = DiscreteCallback(
        (u, t, integrator) -> t == t_resize, affect!;
        save_positions = (false, false)
    )
    return solve(
        ODEProblem(f_resize!, copy(u0_resize), tspan_resize), alg;
        callback = cb, tstops = [t_resize], kwargs...
    )
end

function max_post_resize_error(sol, n_new)
    err = 0.0
    for (t, u) in zip(sol.t, sol.u)
        e = exact_at(t, n_new)
        length(u) == length(e) ||
            error("length mismatch at t=$t: got $(length(u)), expected $(length(e))")
        err = max(err, maximum(abs, u .- e))
    end
    return err
end

function final_resize_error(sol, n_new)
    e = exact_final(n_new)
    length(sol.u[end]) == length(e) ||
        error("final length mismatch: got $(length(sol.u[end])), expected $(length(e))")
    return norm(sol.u[end] - e, Inf)
end

const FIXED_ABM_TYPES = Union{AB3, AB4, AB5, ABM32, ABM43, ABM54}

@testset "resize! grow/shrink: $(nameof(typeof(alg)))" for alg in (
        AB3(), AB4(), AB5(), ABM32(), ABM43(), ABM54(),
        VCAB3(), VCAB4(), VCAB5(), VCABM3(), VCABM4(), VCABM5(), VCABM(),
    )
    fixed = alg isa FIXED_ABM_TYPES
    kwargs = fixed ? (adaptive = false, dt = 0.01) :
        (abstol = 1.0e-10, reltol = 1.0e-10)
    # Coupled exp(At) endpoint residuals with callback restart are ≤3e-7 (AB3) /
    # ≤4e-8 (VCAB*); without restart they jump to ~1e-4–1e-3.
    atol = 1.0e-6

    @testset "grow" begin
        sol = solve_with_resize(alg, 4; kwargs...)
        @test SciMLBase.successful_retcode(sol)
        @test length(sol.u[end]) == 4
        @test final_resize_error(sol, 4) < atol
        @test max_post_resize_error(sol, 4) < 1.0e-5
    end

    @testset "shrink" begin
        sol = solve_with_resize(alg, 2; kwargs...)
        @test SciMLBase.successful_retcode(sol)
        @test length(sol.u[end]) == 2
        @test final_resize_error(sol, 2) < atol
        @test max_post_resize_error(sol, 2) < 1.0e-5
    end
end

@testset "resize! without restart is incorrect on coupled system" begin
    # VCAB3 grow / VCAB5 shrink leave stale ϕ history when u_modified is cleared.
    sol = solve_with_resize(VCAB3(), 4; modified = false, abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test final_resize_error(sol, 4) > 1.0e-5
    sol = solve_with_resize(VCAB5(), 2; modified = false, abstol = 1.0e-10, reltol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol)
    @test final_resize_error(sol, 2) > 1.0e-5
end

@testset "deleteat!: $(nameof(typeof(alg)))" for alg in (
        VCAB3(), VCAB5(), VCABM4(), VCABM(),
    )
    keep = Ref([1, 2, 3])
    A3 = A_resize[1:3, 1:3]
    function f_del!(du, u, p, t)
        mul!(du, view(A3, keep[], keep[]), u)
        return nothing
    end
    function affect!(integrator)
        deleteat!(integrator, 2)
        keep[] = [1, 3]
        return nothing
    end
    cb = DiscreteCallback(
        (u, t, integrator) -> t == t_resize, affect!;
        save_positions = (false, false)
    )
    sol = solve(
        ODEProblem(f_del!, copy(u0_resize), tspan_resize), alg;
        callback = cb, tstops = [t_resize], abstol = 1.0e-10, reltol = 1.0e-10
    )
    u1 = exp(A3 * t_resize) * u0_resize
    ex = exp(A3[[1, 3], [1, 3]] * (tspan_resize[2] - t_resize)) * u1[[1, 3]]
    @test SciMLBase.successful_retcode(sol)
    @test length(sol.u[end]) == 2
    @test norm(sol.u[end] - ex, Inf) < 1.0e-7
end
