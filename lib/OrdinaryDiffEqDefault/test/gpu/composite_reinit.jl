using OrdinaryDiffEqDefault, OrdinaryDiffEqCore, ADTypes, CUDA, Test
import SciMLBase
CUDA.allowscalar(false)

# u = 0 before the step and 1 - exp(-1000 (t - 0.5)) after it, in every component
step_rhs!(du, u, p, t) = (du .= 1000 .* ((t >= 0.5 ? 1.0 : 0.0) .- u); nothing)
step_exact(t) = t <= 0.5 ? 0.0 : -expm1(-1000 * (t - 0.5))
ts = [0.25, 0.5, 0.501, 1.5]

function nlsolver_states(cache::OrdinaryDiffEqCore.DefaultCache)
    subs = (getfield(cache, Symbol(:cache, k)) for k in 1:6 if isdefined(cache, Symbol(:cache, k)))
    return [(s.nlsolver.prev_θ, s.nlsolver.cache.firstcall) for s in subs if hasfield(typeof(s), :nlsolver)]
end

# With reltol < 1e-6 the default algorithm's stiff member is Rodas5P for 8 states and FBDF for 60.
stiff_members = (8 => 4, 60 => 5)
@testset "GPU reinit!: $n states, $(nameof(typeof(ad)))" for (n, stiff) in stiff_members,
        ad in (AutoForwardDiff(), AutoFiniteDiff())

    prob = ODEProblem(step_rhs!, CUDA.zeros(Float64, n), (0.0, 1.5))
    alg = DefaultODEAlgorithm(autodiff = ad)
    kw = (; abstol = 1.0e-8, reltol = 1.0e-8, tstops = [0.5])
    fresh = solve(prob, alg; kw...)
    @test SciMLBase.successful_retcode(fresh)
    @test stiff in fresh.alg_choice
    integ = init(prob, alg; kw...)
    for run in 1:3
        if run > 1
            stiff == 5 && @test any(!=((1.0, true)), nlsolver_states(integ.cache))
            reinit!(integ)
            stiff == 5 && @test all(==((1.0, true)), nlsolver_states(integ.cache))
        end
        solve!(integ)
        @test SciMLBase.successful_retcode(integ.sol)
        @test integ.u isa CuArray
        @test stiff in integ.sol.alg_choice
        @test integ.sol.t == fresh.t
        @test integ.sol.alg_choice == fresh.alg_choice
        @test all(t -> all(x -> isapprox(x, step_exact(t); atol = 1.0e-6, rtol = 1.0e-5), Array(integ.sol(t))), ts)
    end
end
