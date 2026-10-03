using OrdinaryDiffEqPDIRK
using LinearAlgebra
using SciMLBase
using Test

# Coupled linear system: off-diagonal coupling so a wrong stage-buffer size
# cannot accidentally agree with the exact solution on a diagonal problem.
const A = [
    -1.0 0.4 0.2
    0.3 -1.5 0.5
    0.1 0.6 -1.2
]
const U0 = [1.0, -0.5]
const UNEW = 0.3
const TG = 0.5
const TS = 1.0
const TF = 2.0

f!(du, u, p, t) = (mul!(du, view(A, 1:length(u), 1:length(u)), u); nothing)

function exact(t)
    u = exp(A[1:2, 1:2] * min(t, TG)) * U0
    t <= TG && return u
    u3 = exp(A * (min(t, TS) - TG)) * vcat(u, UNEW)
    t <= TS && return u3
    return exp(A[1:2, 1:2] * (t - TS)) * u3[1:2]
end

function resize_callback!()
    grew = Ref(false)
    shrank = Ref(false)
    function condition(u, t, integrator)
        return (!grew[] && t >= TG) || (!shrank[] && t >= TS)
    end
    function affect!(integrator)
        n = grew[] ? 2 : 3
        resize!(integrator, n)
        if n == 3
            integrator.u[3] = UNEW
            # SciML/OrdinaryDiffEq.jl#4722: new `uprev` slots are uninitialized after
            # `resize!`; set them so the solve is deterministic.
            integrator.uprev[3] = UNEW
        end
        grew[] ? (shrank[] = true) : (grew[] = true)
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

const ALG = PDIRK44(threading = false)

@testset "PDIRK44 resize! grow then shrink" begin
    sol = solve(
        ODEProblem(f!, copy(U0), (0.0, TF)), ALG;
        dt = 1.0e-3, adaptive = false, callback = resize_callback!(),
        tstops = [TG, TS]
    )
    @test SciMLBase.successful_retcode(sol)
    @test length(sol.u[end]) == 2
    @test sol(0.75) ≈ exact(0.75) rtol = 1.0e-8 atol = 1.0e-10
    @test sol.u[end] ≈ exact(TF) rtol = 1.0e-8 atol = 1.0e-10
end

# `deleteat!` reaches stage buffers only through `resize_non_user_cache!`
# (default: resize to `length(u)`). Stage values are overwritten each step.
@testset "PDIRK44 deleteat!" begin
    keep = Ref([1, 2, 3])
    A3 = A
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
        (u, t, integrator) -> t == TG, affect!;
        save_positions = (false, false)
    )
    u0 = [1.0, 2.0, 3.0]
    sol = solve(
        ODEProblem(f_del!, copy(u0), (0.0, TF)), ALG;
        dt = 1.0e-3, adaptive = false, callback = cb, tstops = [TG]
    )
    u1 = exp(A3 * TG) * u0
    ex = exp(A3[[1, 3], [1, 3]] * (TF - TG)) * u1[[1, 3]]
    @test SciMLBase.successful_retcode(sol)
    @test length(sol.u[end]) == 2
    @test sol.u[end] ≈ ex rtol = 1.0e-8 atol = 1.0e-10
end
