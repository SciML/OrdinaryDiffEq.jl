using OrdinaryDiffEqBDF, OrdinaryDiffEqCore, Test
using SciMLBase: DAEProblem, successful_retcode, init, step!, solve!
using ADTypes: AutoForwardDiff, AutoFiniteDiff

function dae_lin!(res, du, u, p, t)
    @. res = du - u
    return nothing
end
dae_lin(du, u, p, t) = @. du - u

const prob_iip = DAEProblem(dae_lin!, [1.0, 1.0], [1.0, 1.0], (0.0, 1.0))
const prob_oop = DAEProblem{false}(dae_lin, [1.0, 1.0], [1.0, 1.0], (0.0, 1.0))
const dae_bdf_algs = (DImplicitEuler(), DABDF2(), DFBDF(), DNordsieckBDF())

@testset "DAE BDF uses COEFFICIENT_MULTISTEP (issue #4113)" begin
    for (prob, label) in ((prob_iip, "iip"), (prob_oop, "oop"))
        @testset "$label" begin
            for alg in dae_bdf_algs
                integrator = init(
                    prob, alg; dt = 0.1, adaptive = false, save_everystep = false
                )
                step!(integrator)
                step!(integrator)
                nlsolver = integrator.cache.nlsolver
                @test nlsolver.method === OrdinaryDiffEqCore.COEFFICIENT_MULTISTEP
                @test nlsolver.cache.W_γdt ≈
                    nlsolver.α / (nlsolver.γ * integrator.dt)
            end
        end
    end
end

@testset "DFBDF startup retries use the accepted state" begin
    function residual!(r, du, u, p, t)
        exponent, shift = p
        factor = exponent == 0 ? one(du[1]) : 1 + (du[1] + shift)^exponent + u[2]^2
        r[1] = (du[1] + u[1]) * factor
        r[2] = u[2] - u[1]^2
        return nothing
    end
    function residual(du, u, p, t)
        r = similar(u, promote_type(eltype(du), eltype(u)))
        residual!(r, du, u, p, t)
        return r
    end
    for p in ((0, 0.0), (2, 0.0), (4, 0.0), (2, 0.5)), iip in (true, false),
            autodiff in (AutoForwardDiff(), AutoFiniteDiff()), initial_dt in (nothing, 1.0e-8)
        @testset "p=$p iip=$iip autodiff=$autodiff dt=$initial_dt" begin
            prob = DAEProblem(
                iip ? residual! : residual, [-1.0, -2.0], [1.0, 1.0], (0.0, 1.0), p;
                differential_vars = [true, false]
            )
            dt_kwargs = initial_dt === nothing ? NamedTuple() : (; dt = initial_dt)
            integrator = init(
                prob, DFBDF(; autodiff); abstol = 1.0e-7, reltol = 1.0e-7, dt_kwargs...
            )
            step!(integrator)
            @test integrator.t > 0
            @test initial_dt !== nothing || integrator.sol.stats.nreject > 0
            @test integrator.sol.u[1] == [1.0, 1.0]
            sol = solve!(integrator)
            @test successful_retcode(sol)
            @test sol.t[end] == 1.0
            @test sol.u[1] == [1.0, 1.0]
            @test maximum(
                maximum(abs.(u .- [exp(-t), exp(-2t)])) for (u, t) in zip(sol.u, sol.t)
            ) < 2.0e-6
            @test maximum(abs(u[2] - u[1]^2) for u in sol.u) < 1.0e-6
            @test prob.u0 == [1.0, 1.0]
            @test prob.du0 == [-1.0, -2.0]
        end
    end
end
