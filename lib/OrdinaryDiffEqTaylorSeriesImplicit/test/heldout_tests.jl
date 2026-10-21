# Problems not used to develop the solver, checked for properties that hold for any problem:
#   * every solve succeeds with an endpoint error within a few tolerance units;
#   * tightening the tolerance does not increase the error (beyond roundoff);
#   * repeated solves are bit-identical;
#   * fixed-step solves converge with the classical order.
using Test, LinearAlgebra
using OrdinaryDiffEqTaylorSeriesImplicit, OrdinaryDiffEqRosenbrock
import SciMLBase

module HeldOut
using LinearAlgebra, SciMLBase
# Prothero-Robinson, u = sin t
pr!(du, u, p, t) = (du[1] = p * (u[1] - sin(t)) + cos(t); nothing)
# Kaps, u = (e^{-2t}, e^{-t}) for any ε
kaps!(du, u, p, t) = (du[1] = -(1 / p + 2) * u[1] + u[2]^2 / p; du[2] = u[1] - u[2] - u[2]^2; nothing)
# Stiff linear system with a complex pair, u = exp(At) u0
const A = let V = [1.0 0.3 0.1 0.0; 0.2 1.0 0.0 0.4; 0.0 0.1 1.0 0.2; 0.3 0.0 0.2 1.0]
    V * [-1.0 0 0 0; 0 -2e3 300.0 0; 0 -300.0 -2e3 0; 0 0 0 -1e5] / V
end
lin!(du, u, p, t) = (mul!(du, A, u); nothing)
# Curtiss-Hirschfelder
ch!(du, u, p, t) = (du[1] = -50 * (u[1] - cos(t)); nothing)
ch_exact(t) = 2500 / 2501 * cos(t) + 50 / 2501 * sin(t) - 2500 / 2501 * exp(-50t)
# Reaction-diffusion u_t = D u_xx + 5u(1-u)(u-a) with zero Dirichlet data
function rd!(du, u, p, t)
    N = length(u); D, a = p; h2 = (1 / (N + 1))^2
    @inbounds for i in 1:N
        l = i == 1 ? 0.0 : u[i - 1]; r = i == N ? 0.0 : u[i + 1]
        du[i] = D * (l - 2u[i] + r) / h2 + 5 * u[i] * (1 - u[i]) * (u[i] - a)
    end
end
rd(N, T) = ODEProblem(rd!, [sin(pi * i / (N + 1))^2 for i in 1:N], (0.0, T), (0.5, 0.3))
# (name, problem, exact endpoint or nothing for a Rodas5P reference)
problems() = [
    ("Prothero-Robinson", ODEProblem(pr!, [0.0], (0.0, 10.0), -1e4), [sin(10.0)]),
    ("Kaps eps=1e-6", ODEProblem(kaps!, [1.0, 1.0], (0.0, 5.0), 1e-6), [exp(-10.0), exp(-5.0)]),
    ("stiff linear, complex pair", ODEProblem(lin!, [1.0, 0.5, -0.3, 0.8], (0.0, 3.0)), exp(3A) * [1.0, 0.5, -0.3, 0.8]),
    ("Curtiss-Hirschfelder", ODEProblem(ch!, [0.0], (0.0, 3.0)), [ch_exact(3.0)]),
    ("reaction-diffusion n=20", rd(20, 2.0), nothing),
    ("reaction-diffusion n=40", rd(40, 2.0), nothing),
    ("Kaps backward in time", ODEProblem(kaps!, [exp(-2.0), exp(-1.0)], (1.0, 0.0), 1.0), [1.0, 1.0]),
]
end

@testset "Implicit Taylor held-out problems" begin
    reference(prob, ex) = ex !== nothing ? ex :
        solve(prob, Rodas5P(); abstol = 1e-14, reltol = 1e-13, save_everystep = false).u[end]
    scaled(u, ref, at, rt) = sqrt(sum(abs2, (u .- ref) ./ (at .+ rt .* abs.(ref))) / length(u))
    run(prob, alg, rt) = solve(prob, alg; abstol = rt / 100, reltol = rt, save_everystep = false, maxiters = 10^6)
    methods = (("GaussA4", TaylorGaussA4), ("RadauLA3", TaylorRadauLA3), ("RadauLA5", TaylorRadauLA5),
        ("LobattoLA4", TaylorLobattoLA4), ("LobattoLA6", TaylorLobattoLA6))

    @testset "$name" for (name, prob, ex) in HeldOut.problems()
        ref = reference(prob, ex)
        floor = 1e-10 * max(1.0, maximum(abs, ref))      # roundoff level of the endpoint
        @testset "$mname" for (mname, M) in methods
            errs = Float64[]
            for rt in (1e-4, 1e-7)
                s = run(prob, M(), rt)
                @test SciMLBase.successful_retcode(s)
                @test scaled(s.u[end], ref, rt / 100, rt) < 5
                push!(errs, norm(s.u[end] .- ref, Inf))
                # bit-identical repeat
                s2 = run(prob, M(), rt)
                @test s2.u[end] == s.u[end]
                @test (s2.stats.naccept, s2.stats.nw, s2.stats.nnonliniter) == (s.stats.naccept, s.stats.nw, s.stats.nnonliniter)
            end
            # Not 'decreases': a loose-tolerance solve may already be far more accurate than asked.
            @test errs[2] <= max(errs[1], floor)
        end
    end

    @testset "Fixed-step order on Kaps (eps = 1)" begin
        prob = ODEProblem(HeldOut.kaps!, [1.0, 1.0], (0.0, 1.0), 1.0)
        exact = [exp(-2.0), exp(-1.0)]
        for (M, order, dts) in ((TaylorGaussA4, 4, 2.0 .^ (-3:-1:-6)), (TaylorRadauLA3, 3, 2.0 .^ (-3:-1:-6)),
                (TaylorRadauLA5, 5, 2.0 .^ (-2:-1:-5)), (TaylorLobattoLA4, 4, 2.0 .^ (-3:-1:-6)),
                (TaylorLobattoLA6, 6, 2.0 .^ (-1:-1:-4)))
            errs = [norm(solve(prob, M(); dt, adaptive = false, save_everystep = false).u[end] .- exact, Inf) for dt in dts]
            slopes = diff(log2.(errs)) ./ diff(log2.(dts))
            @test slopes[end] ≈ order atol = 0.3
        end
    end
end
