using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, ForwardDiff, SciMLBase, Test
using OrdinaryDiffEqCore: PIController, zero_func_struct

# u' = p for t < 1/2 and u' = 2p after, u(0) = 0: u(1) = 3p/2, du(1)/dp = 3/2.
kink!(du, u, p, t) = (du[1] = (t < 0.5 ? 1 : 2) * p[1]; nothing)

const ZERO_FUNC_CALLS = Ref(0)
count_zero_func(integ) = (integ isa zero_func_struct && (ZERO_FUNC_CALLS[] += 1); nothing)
scond(u, t, integ) = (count_zero_func(integ); t - 0.5)
function vcond(out, u, t, integ)
    count_zero_func(integ)
    out[1] = t - 0.5
    out[2] = u[1] + 1
    return nothing
end

function u_end(p, spec, cb)
    prob = ODEProblem{true, spec}(kink!, [zero(p)], (0.0, 1.0), [p])
    sol = solve(
        prob, Tsit5(); callback = cb, dt = 0.8, reltol = 1.0e-10, abstol = 1.0e-12,
        controller = PIController(Tsit5(); discontinuity_detection = true)
    )
    return sol.u[end][1]
end

@testset "discontinuity detection with ForwardDiff duals" begin
    p = 1.3
    callbacks = (
        ("scalar", ContinuousCallback(scond, integ -> nothing)),
        ("vector", VectorContinuousCallback(vcond, (integ, i) -> nothing, 2)),
    )
    @testset "$name, $spec" for (name, cb) in callbacks,
            spec in (SciMLBase.FullSpecialize, SciMLBase.AutoSpecialize)
        ZERO_FUNC_CALLS[] = 0
        @test u_end(p, spec, cb) ≈ 1.5p rtol = 1.0e-6
        @test ZERO_FUNC_CALLS[] > 0
        ZERO_FUNC_CALLS[] = 0
        @test ForwardDiff.derivative(q -> u_end(q, spec, cb), p) ≈ 1.5 rtol = 1.0e-6
        @test ZERO_FUNC_CALLS[] > 0
    end
end
