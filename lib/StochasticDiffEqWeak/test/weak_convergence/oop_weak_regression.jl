"""
Regression tests for previously-broken out-of-place weak SDE paths.

These exercise the constant-cache (out-of-place) `perform_step!` code paths for
PL1WM (Itô) and RS1/RS2 (Stratonovich) with multi-dimensional diagonal and
non-diagonal noise. Every one of these used to error before the fixes:

  * PL1WM / RS `H12`, `Yp`, `Ym` were built with `Vector{typeof(uprev)}[...]`,
    which wraps the element type one level too deep -> `convert` `MethodError`.
  * RS `H13`/`H14` non-diagonal updates dropped the per-channel `[k]` index and
    seeded the stage vectors with `zero(typeof(uprev))` (undefined for arrays).
  * RS diagonal cross-terms used `g2[l][l]` (scalar) instead of `g2[l]`.
  * RS `H22`/`H23` aliased `uprev` (`[uprev for k in 1:m]`), so the in-place
    accumulation corrupted `uprev` itself.

The mean of a linear SDE with diagonal drift `du_i = a_i u_i dt + noise` is known
in closed form, so we check weak accuracy against it:
  * Itô (PL1WM):          E[u_i(T)] = u0_i * exp(a_i T)
  * Stratonovich (RS):    E[u_i(T)] = u0_i * exp((a_i + 1/2 * sum_k b_{i,k}^2) T)
"""

using StochasticDiffEqWeak
using Test
using Random
import Statistics

const T_END = 1.0
const U0 = [1.0, 2.0]
const AVEC = [-0.5, -0.3]                 # diagonal drift coefficients
# diffusion coefficients b[i, k] = d(noise channel k of component i)/du_i
const BDIAG = [0.3 0.0; 0.0 0.25]         # diagonal noise (channel k drives only component k)
const BND = [0.3 0.1; 0.15 0.25]          # non-diagonal noise matrix g(u)[i,k] = b[i,k]*u[i]

ito_mean() = [U0[i] * exp(AVEC[i] * T_END) for i in 1:2]
function strat_mean(b)
    return [U0[i] * exp((AVEC[i] + 0.5 * sum(b[i, k]^2 for k in 1:2)) * T_END) for i in 1:2]
end

# out-of-place problem definitions
fd(u, p, t) = [AVEC[1] * u[1], AVEC[2] * u[2]]
gd(u, p, t) = [BDIAG[1, 1] * u[1], BDIAG[2, 2] * u[2]]
gn(u, p, t) = [BND[1, 1] * u[1] BND[1, 2] * u[1]; BND[2, 1] * u[2] BND[2, 2] * u[2]]

prob_diag = SDEProblem(fd, gd, copy(U0), (0.0, T_END))
prob_nd = SDEProblem(fd, gn, copy(U0), (0.0, T_END), noise_rate_prototype = zeros(2, 2))

function weak_mean(prob, alg, dt, N)
    ep = EnsembleProblem(prob)
    sol = solve(
        ep, alg, EnsembleThreads(); trajectories = N, dt = dt,
        save_everystep = false, save_start = false, adaptive = false
    )
    us = [sol.u[i].u[end] for i in 1:N]
    return [Statistics.mean(u[j] for u in us) for j in 1:2]
end

const N = 60_000
const DT = 1 // 50
const TOL = 0.02

@testset "out-of-place weak SDE regression" begin
    @testset "PL1WM (Itô) out-of-place non-diagonal" begin
        Random.seed!(20)
        m = weak_mean(prob_nd, PL1WM(), DT, N)
        @test all(isfinite, m)
        @test maximum(abs, m .- ito_mean()) < TOL
    end

    for alg in (RS1(), RS2())
        nm = string(nameof(typeof(alg)))
        @testset "$nm (Stratonovich) out-of-place diagonal" begin
            u0_before = copy(prob_diag.u0)
            Random.seed!(20)
            m = weak_mean(prob_diag, alg, DT, N)
            @test all(isfinite, m)
            @test maximum(abs, m .- strat_mean(BDIAG)) < TOL
            @test prob_diag.u0 == u0_before   # an out-of-place solve must not mutate u0
        end
        @testset "$nm (Stratonovich) out-of-place non-diagonal" begin
            u0_before = copy(prob_nd.u0)
            Random.seed!(20)
            m = weak_mean(prob_nd, alg, DT, N)
            @test all(isfinite, m)
            @test maximum(abs, m .- strat_mean(BND)) < TOL
            @test prob_nd.u0 == u0_before
        end
    end
end
