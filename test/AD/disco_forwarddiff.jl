using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, ForwardDiff, SciMLBase, Test
using OrdinaryDiffEqCore: PIController, zero_func_struct

# u' = p for t < 1/2 and u' = 2p after, u(0) = 0: u(1) = 3p/2.
kink!(du, u, p, t) = (du[1] = (t < 0.5 ? 1 : 2) * p[1]; nothing)
# u' = p while u < 1/2 and u' = 2p after, u(0) = 0: the crossing is at t = 1/(2p), u(1) = 2p - 1/2.
state_kink!(du, u, p, t) = (du[1] = (u[1] < 0.5 ? 1 : 2) * p[1]; nothing)

const ZERO_FUNC_CALLS = Ref(0)
count_zero_func(integ) = (integ isa zero_func_struct && (ZERO_FUNC_CALLS[] += 1); nothing)
scalar_condition(u, t, integ) = (count_zero_func(integ); t - 0.5)
function vcond(out, u, t, integ)
    count_zero_func(integ)
    out[1] = t - 0.5
    out[2] = u[1] + 1
    return nothing
end
state_condition(u, t, integ) = (count_zero_func(integ); u[1] - 0.5)

function u_end(f!, p, spec, cb)
    prob = ODEProblem{true, spec}(f!, [zero(p)], (0.0, 1.0), [p])
    sol = solve(
        prob, Tsit5(); callback = cb, dt = 0.8, reltol = 1.0e-10, abstol = 1.0e-12,
        controller = PIController(Tsit5(); discontinuity_detection = true)
    )
    return sol.u[end][1]
end

# Discontinuity detection does not differentiate through the located root: a Dual state
# reaching it throws converting the Dual to Float64 instead of returning a derivative.
@testset "discontinuity detection: Float64 solves and ForwardDiff duals" begin
    p = 1.3
    cases = (
        ("scalar", kink!, ContinuousCallback(scalar_condition, integ -> nothing), 1.5p),
        ("vector", kink!, VectorContinuousCallback(vcond, (integ, i) -> nothing, 2), 1.5p),
        ("state", state_kink!, ContinuousCallback(state_condition, integ -> nothing), 2p - 0.5),
    )
    @testset "$name, $spec" for (name, f!, cb, expected) in cases,
            spec in (SciMLBase.FullSpecialize, SciMLBase.AutoSpecialize)
        ZERO_FUNC_CALLS[] = 0
        @test u_end(f!, p, spec, cb) ≈ expected rtol = 1.0e-6
        @test ZERO_FUNC_CALLS[] > 0
        err = try
            ForwardDiff.derivative(q -> u_end(f!, q, spec, cb), p)
        catch e
            e
        end
        @test err isa MethodError && err.f === Float64 && only(err.args) isa ForwardDiff.Dual
    end
end
