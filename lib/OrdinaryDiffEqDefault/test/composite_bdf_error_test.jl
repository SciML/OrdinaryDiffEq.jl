using OrdinaryDiffEqDefault, OrdinaryDiffEqBDF, OrdinaryDiffEqCore, OrdinaryDiffEqTsit5,
    OrdinaryDiffEqVerner, SciMLBase, Test

# Fast states relax to a forcing that jumps at `tjump`, and slow states integrate them, so a jump
# taken at the wrong time leaves a permanent error in the slow states:
#   x' = λ (H(t - tjump) - x), x(0) = 1
#   y' = x,                    y(0) = 0
const NPAIRS = 30
const tjump = 1.0
const tend = 2.0
rate(i) = 1000.0 * (1 + (i - 1) / NPAIRS)

function jump_rhs!(du, u, p, t)
    g = t >= tjump ? 1.0 : 0.0
    for i in 1:NPAIRS
        du[i] = rate(i) * (g - u[i])
        du[NPAIRS + i] = u[i]
    end
    return nothing
end

exact_x(λ, t) = t < tjump ? exp(-λ * t) :
    1 + (exp(-λ * tjump) - 1) * exp(-λ * (t - tjump))
function exact_y(λ, t)
    t < tjump && return (1 - exp(-λ * t)) / λ
    xj = exp(-λ * tjump)
    return (1 - xj) / λ + (t - tjump) + (xj - 1) * (1 - exp(-λ * (t - tjump))) / λ
end

prob = ODEProblem(jump_rhs!, vcat(ones(NPAIRS), zeros(NPAIRS)), (0.0, tend))

@testset "A BDF member of a composite honours its error test: $name" for (name, alg) in (
            "DefaultODEAlgorithm" => DefaultODEAlgorithm(),
            "AutoTsit5(FBDF)" => AutoTsit5(FBDF()),
            "AutoVern7(FBDF)" => AutoVern7(FBDF()),
            "AutoVern9(FBDF)" => AutoVern9(FBDF()),
        ), tstops in (Float64[], [tjump]), (abstol, reltol) in ((1.0e-6, 1.0e-3), (1.0e-8, 1.0e-6))

    integ = init(prob, alg; abstol, reltol, tstops)
    accepted_over_tolerance = 0
    while integ.t < tend
        step!(integ)
        accepted_over_tolerance += OrdinaryDiffEqCore.get_EEst(integ) > 1
    end
    solve!(integ)
    sol = integ.sol

    @test sol.retcode == ReturnCode.Success
    @test any(>(1), sol.alg_choice)
    @test accepted_over_tolerance == 0
    for i in (1, NPAIRS)
        λ = rate(i)
        @test sol.u[end][NPAIRS + i] ≈ exact_y(λ, tend) rtol = 10 * reltol
        @test sol(0.99)[i] ≈ exact_x(λ, 0.99) atol = 100 * abstol
    end
end
