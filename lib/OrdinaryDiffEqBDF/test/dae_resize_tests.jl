using OrdinaryDiffEqBDF, LinearAlgebra, ADTypes, Test

# Global components 1:3 are differential, 4:5 algebraic. The active components, in state
# order, are ACTIVE[]. Differential rows: x_k' = M[k, active] ⋅ u. Algebraic rows:
# y_k = M[k, active] ⋅ u, where M has no algebraic columns in algebraic rows (index 1).
const M = [
    -1.0 0.3 0.1 0.2 -0.1
    0.2 -1.4 0.25 0.1 0.15
    0.15 0.3 -1.2 -0.2 0.1
    0.5 -0.4 0.3 0.0 0.0
    -0.2 0.6 0.4 0.0 0.0
]
const ISDIFF = [true, true, true, false, false]
const ACTIVE = Ref([1, 2, 4])
const TG = 0.5
const TF = 2.0

function residual!(r, du, u, p, t)
    act = ACTIVE[]
    mul!(r, view(M, act, act), u)
    for (i, k) in enumerate(act)
        r[i] = (ISDIFF[k] ? du[i] : u[i]) - r[i]
    end
    return nothing
end

# Make `u` consistent for the components `act`: the algebraic entries are solved from the
# (index-1) constraints given the differential entries, and `du` is the exact derivative.
function consistent!(du, u, act)
    d = findall(k -> ISDIFF[k], act)
    a = findall(k -> !ISDIFF[k], act)
    u[a] = M[act[a], act[d]] * u[d]
    du[d] = M[act[d], act] * u
    du[a] = M[act[a], act[d]] * du[d]
    return nothing
end

function exact(u0, act, t)
    d = findall(k -> ISDIFF[k], act)
    a = findall(k -> !ISDIFF[k], act)
    C = M[act[a], act[d]]
    x = exp((M[act[d], act[d]] + M[act[d], act[a]] * C) * t) * u0[d]
    u = similar(u0)
    u[d] = x
    u[a] = C * x
    return u
end

@testset "DFBDF resize!/deleteat!/addat! mid-solve" begin
    # `addat!` on arrays appends, so the new component is last.
    cases = (
        ("grow", [1, 2, 4, 5], integ -> resize!(integ, 4)),
        ("shrink", [1, 2], integ -> resize!(integ, 2)),
        ("deleteat!", [1, 4], integ -> deleteat!(integ, [2])),
        ("addat!", [1, 2, 4, 3], integ -> (addat!(integ, [2]); integ.u[4] = 0.3)),
    )
    # A ForwardDiff chunk larger than the shrunk state is a separate limitation (#4803).
    algs = (
        "ForwardDiff" => DFBDF(autodiff = AutoForwardDiff(chunksize = 1)),
        "FiniteDiff" => DFBDF(autodiff = AutoFiniteDiff()),
    )
    for (adname, alg) in algs, (name, after, edit!) in cases
        @testset "$name $adname" begin
            before = [1, 2, 4]
            ACTIVE[] = before
            u0 = [1.0, -0.5, 0.0]
            du0 = zeros(3)
            consistent!(du0, u0, before)
            prob = DAEProblem(
                residual!, du0, u0, (0.0, TF);
                differential_vars = ISDIFF[before]
            )
            integ = init(prob, alg; tstops = [TG], reltol = 1.0e-10, abstol = 1.0e-12)
            while integ.t < TG
                step!(integ)
            end
            @test integ.u ≈ exact(u0, before, TG) rtol = 1.0e-7

            edit!(integ)
            ACTIVE[] = after
            @test length(integ.u) == length(integ.du) == length(after)
            # New algebraic components, and old ones whose constraint changed, are set
            # from the constraints so the edited state is consistent.
            consistent!(integ.du, integ.u, after)
            u_edit = copy(integ.u)
            du_edit = copy(integ.du)
            r = similar(u_edit)
            residual!(r, du_edit, u_edit, nothing, TG)
            @test norm(r) < 1.0e-12

            solve!(integ)
            @test integ.sol.retcode == ReturnCode.Success

            fresh = solve(
                DAEProblem(
                    residual!, du_edit, u_edit, (TG, TF);
                    differential_vars = ISDIFF[after]
                ), alg; reltol = 1.0e-10, abstol = 1.0e-12
            )
            @test fresh.retcode == ReturnCode.Success
            @test integ.u ≈ fresh.u[end] rtol = 1.0e-6
            @test integ.u ≈ exact(u_edit, after, TF - TG) rtol = 1.0e-6
        end
    end
end
