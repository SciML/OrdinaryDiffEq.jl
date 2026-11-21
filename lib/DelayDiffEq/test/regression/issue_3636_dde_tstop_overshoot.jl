using Test, DelayDiffEq, DiffEqCallbacks, OrdinaryDiffEqTsit5
using SciMLBase: ReturnCode

# Regression test for SciML/OrdinaryDiffEq.jl#3636.
#
# With a neutral DDE whose constant lag is an exact integer multiple of dt,
# `tprev + clamped_dt` can round one ULP above the next tstop. The ODE path
# guards against this via `next_step_tstop` / `tstop_target` snapping in
# `fixed_t_for_tstop_error!`; before this fix `DDEIntegrator` did not carry
# those fields, so the fallback returned `ttmp` unchanged and `handle_tstop!`
# raised "Integrator stepped past tstops but the algorithm was dtchangeable".
#
# Three-state partition for G(s) = 1/(0.85s+1) * exp(-0.14s) driven by u(t)=1:
#   dx = A*x + d(t-τ),  dY = C1*x,  dD = 1
# The neutral lag pulls d(t-τ) as the Val{1} derivative of history at idx 3.

const A_ = -1.0 / 0.85
const C1_ = 1.0 / 0.85

function dde_3636!(du, u, h, p, t)
    τ = p[1]
    x = u[1]
    d_lag = h(p, t - τ, Val{1}; idxs = 3)
    du[1] = A_ * x + d_lag
    du[2] = C1_ * x
    du[3] = 1.0
    return nothing
end

h_3636(p, t, ::Type{Val{1}}; idxs = 0) = 0.0

@testset "issue #3636: DDE tstop overshoot by one ULP" begin
    τ = 0.14
    # The original failure mode is dt = 0.014 (τ / dt = 10). Sweep neighbouring
    # values too so that we exercise both the snap-back path and the regular
    # path.
    for dt in (0.011, 0.012, 0.013, 0.0125, 0.0117, 0.014, 0.015, 0.016, 0.02, 0.028)
        tgrid = 0:dt:5.0
        prob = DDEProblem{true}(
            dde_3636!, [0.0, 0.0, 0.0], h_3636, (0.0, 5.0), (τ,);
            constant_lags = [τ], neutral = true
        )
        sv = SavedValues(Float64, Tuple{Vector{Float64}, Vector{Float64}})
        tmpy = [0.0]
        saver(u, t, integ) = (u[1:1], copy(tmpy))
        cb = SavingCallback(saver, sv; saveat = tgrid)

        sol = solve(
            prob, MethodOfSteps(Tsit5());
            tstops = tgrid, saveat = tgrid,
            abstol = 1.0e-6, reltol = 1.0e-6, force_dtmin = true,
            callback = cb
        )
        @test sol.retcode == ReturnCode.Success
        @test !isempty(sv.t)
    end
end
