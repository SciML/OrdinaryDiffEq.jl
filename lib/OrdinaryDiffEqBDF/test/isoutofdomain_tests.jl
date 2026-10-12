using OrdinaryDiffEqBDF, Test
using SciMLBase: DAEProblem, ODEProblem, ReturnCode

# Regression tests for https://github.com/SciML/OrdinaryDiffEq.jl/issues/4687:
# a single `isoutofdomain` rejection must leave the multistep history consistent,
# so the retry with a smaller `dt` is as accurate as an ordinary step.

const TSPAN = (0.0, 40.0)
const TOL = 1.0e-6
const T_REJECT = 5.0

exact(t) = [exp(-0.5 * t), exp(-0.1 * t)]

f_oop(u, p, t) = [-0.5 * u[1], -0.1 * u[2]]
function f_iip(du, u, p, t)
    du[1] = -0.5 * u[1]
    du[2] = -0.1 * u[2]
    return nothing
end

res_oop(du, u, p, t) = du .- f_oop(u, p, t)
function res_iip(r, du, u, p, t)
    r[1] = du[1] + 0.5 * u[1]
    r[2] = du[2] + 0.1 * u[2]
    return nothing
end

ode_problems() = (
    ("in-place", ODEProblem(f_iip, [1.0, 1.0], TSPAN)),
    ("out-of-place", ODEProblem{false}(f_oop, [1.0, 1.0], TSPAN)),
)

dae_problems() = (
    (
        "in-place",
        DAEProblem(res_iip, [-0.5, -0.1], [1.0, 1.0], TSPAN; differential_vars = [true, true]),
    ),
    (
        "out-of-place",
        DAEProblem{false}(
            res_oop, [-0.5, -0.1], [1.0, 1.0], TSPAN; differential_vars = [true, true]
        ),
    ),
)

# Reject exactly one proposed step past `T_REJECT`, as a positivity check would.
function solve_rejecting_once(prob, alg; reject_once)
    nrejected = Ref(0)
    function isoutofdomain(u, p, t)
        if reject_once && t > T_REJECT && nrejected[] == 0
            nrejected[] += 1
            return true
        end
        return false
    end
    sol = solve(prob, alg; isoutofdomain, abstol = TOL, reltol = TOL)
    err = maximum(abs.(sol.u[end] .- exact(sol.t[end])))
    return sol, err, nrejected[]
end

function check_single_domain_rejection(prob, alg)
    sol_ref, err_ref, _ = solve_rejecting_once(prob, alg; reject_once = false)
    sol, err, nrejected = solve_rejecting_once(prob, alg; reject_once = true)
    @test nrejected == 1
    @test sol_ref.retcode == ReturnCode.Success
    @test sol.retcode == ReturnCode.Success
    @test sol.t[end] == TSPAN[2]
    # Before the fix NordsieckBDF lost ~3 orders of magnitude and FBDF aborted with
    # dt below eps. One extra rejection should not change the error materially.
    @test err <= 10 * max(err_ref, 10 * eps())
    @test err <= 100 * TOL
    return nothing
end

@testset "isoutofdomain rejection keeps ODE multistep history consistent (#4687)" begin
    for alg in (NordsieckBDF(), FBDF(), QNDF()), (nm, prob) in ode_problems()
        @testset "$(nameof(typeof(alg))) $nm" begin
            check_single_domain_rejection(prob, alg)
        end
    end
end

@testset "isoutofdomain rejection keeps DAE multistep history consistent (#4687)" begin
    for alg in (DNordsieckBDF(), DFBDF()), (nm, prob) in dae_problems()
        @testset "$(nameof(typeof(alg))) $nm" begin
            check_single_domain_rejection(prob, alg)
        end
    end
end

@testset "isoutofdomain rejection still shrinks dt for multistep methods" begin
    # An out-of-domain step can have a tiny error estimate, so the retry must not be
    # chosen from EEst: every rejection has to reduce dt.
    for alg in (NordsieckBDF(), FBDF(), QNDF())
        @testset "$(nameof(typeof(alg)))" begin
            armed = Ref(false)
            ntries = Ref(0)
            proposed_t = Float64[]
            function isoutofdomain(u, p, t)
                armed[] || return false
                ntries[] += 1
                push!(proposed_t, t)
                return ntries[] <= 3  # reject three consecutive attempts
            end
            prob = ODEProblem(f_iip, [1.0, 1.0], TSPAN)
            integ = init(prob, alg; isoutofdomain, abstol = TOL, reltol = TOL)
            while integ.t <= T_REJECT
                step!(integ)
            end
            t0 = integ.t
            armed[] = true
            step!(integ)
            @test integ.t > t0
            @test length(proposed_t) >= 4
            dts = proposed_t .- t0
            @test all(dts[i + 1] < dts[i] for i in 1:3)
            @test all(>(0), dts)
        end
    end
end
