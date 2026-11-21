using OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqDifferentiation: JVPCache, set_jvp_point!, sync_jvp_point!
using ADTypes, SciMLBase, LinearAlgebra, AllocCheck, Test
using SciMLOperators: update_coefficients!

# `JVPCache`'s product path is non-allocating: fixing the linearization point is
# proven statically with `check_allocs`, and `mul!` is asserted at run time.
# `mul!`'s residual `check_allocs` sites all live inside the DI backend
# `pushforward!` internals — conditional `unaliascopy` sites that never fire for
# the distinct preallocated buffers used here — so those are tracked as broken
# until DifferentiationInterface sheds them.
const N = 24
const dx = 1 / (N + 1)
function allencahn!(du, u, p, t)
    for i in 1:N
        um = i == 1 ? -one(eltype(u)) : u[i - 1]
        up = i == N ? one(eltype(u)) : u[i + 1]
        du[i] = 1.0e-3 * (um - 2u[i] + up) / dx^2 + u[i] - u[i]^3
    end
    return nothing
end
function allencahn_jvp!(Jv, v, u, p, t)
    for i in 1:N
        um = i == 1 ? -one(eltype(v)) : v[i - 1]
        up = i == N ? one(eltype(v)) : v[i + 1]
        Jv[i] = 1.0e-3 * (um - 2v[i] + up) / dx^2 + v[i] - 3u[i]^2 * v[i]
    end
    return nothing
end
u0 = [tanh((i * dx - 0.5) / 0.1) for i in 1:N]
prob = ODEProblem(allencahn!, u0, (0.0, 1.0))
prob_jvp = ODEProblem(ODEFunction(allencahn!; jvp = allencahn_jvp!), u0, (0.0, 1.0))

function alloc_mul_steady(Jv, J, v)
    mul!(Jv, J, v)
    return @allocated mul!(Jv, J, v)
end
function alloc_mul_fixing(Jv, J, v, u, p, t)
    update_coefficients!(J, u, p, t)
    mul!(Jv, J, v)
    update_coefficients!(J, u, p, t)
    return @allocated mul!(Jv, J, v)
end

@testset "JVPCache allocation checks" begin
    @testset "$adname" for (adname, ad) in (
            ("AutoFiniteDiff", AutoFiniteDiff()), ("AutoForwardDiff", AutoForwardDiff()),
        )
        J = JVPCache(prob.f, similar(u0), u0, prob.p, 0.0; autodiff = ad)
        v, Jv = ones(N), zeros(N)
        update_coefficients!(J, u0, prob.p, 0.0)
        mul!(Jv, J, v)
        op = J.jvp_op

        @testset "point-fixing functions are allocation-free" begin
            @test isempty(
                check_allocs(
                    update_coefficients!,
                    (typeof(J), typeof(u0), typeof(prob.p), Float64)
                )
            )
            @test isempty(check_allocs(sync_jvp_point!, (typeof(J),)))
            @test isempty(
                check_allocs(
                    set_jvp_point!,
                    (typeof(op), typeof(u0), typeof(prob.p), Float64)
                )
            )
        end

        # The user-visible guarantee: zero bytes per product, including the
        # product that fixes the point after `update_coefficients!`.
        @test alloc_mul_steady(Jv, J, v) == 0
        @test alloc_mul_fixing(Jv, J, v, u0, prob.p, 0.0) == 0

        allocs = check_allocs(
            op, (typeof(Jv), typeof(v), typeof(u0), typeof(prob.p), Float64)
        )
        @test isempty(allocs) broken = true
        if !isempty(allocs)
            println(
                "AllocCheck found $(length(allocs)) conditional sites in " *
                    "the $adname pushforward! path (inside DifferentiationInterface)"
            )
        end
    end

    # A user-supplied `f.jvp` applies per product and is statically clean end to end.
    @testset "user-supplied f.jvp" begin
        J = JVPCache(prob_jvp.f, similar(u0), u0, prob_jvp.p, 0.0; autodiff = AutoFiniteDiff())
        v, Jv = ones(N), zeros(N)
        update_coefficients!(J, u0, prob_jvp.p, 0.0)
        mul!(Jv, J, v)

        @test alloc_mul_steady(Jv, J, v) == 0
        @test alloc_mul_fixing(Jv, J, v, u0, prob_jvp.p, 0.0) == 0
        @test isempty(
            check_allocs(LinearAlgebra.mul!, (typeof(Jv), typeof(J), typeof(v)))
        )
    end
end
