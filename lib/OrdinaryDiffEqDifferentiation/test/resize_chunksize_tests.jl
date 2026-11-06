using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, SciMLBase, ADTypes, Test
using OrdinaryDiffEqDifferentiation
using LinearAlgebra: diag
import DifferentiationInterface as DI

# The ForwardDiff chunk picked at `init` for 10 states is larger than the shrunk state.
function f!(du, u, p, t)
    n = length(u)
    for i in 1:n
        du[i] = -i * u[i] + u[mod1(i + 1, n)] / 4
    end
    return nothing
end

resize_shrink!(integ) = resize!(integ, 2)
deleteat_shrink!(integ) = deleteat!(integ, 3:10)
cases = (
    (Rodas5P(), resize_shrink!, false),
    (Rodas5P(), deleteat_shrink!, false),
    (Rodas5P(autodiff = AutoForwardDiff(chunksize = 5)), resize_shrink!, false),
    (Rodas5P(autodiff = AutoForwardDiff(chunksize = 5)), deleteat_shrink!, false),
    (TRBDF2(), resize_shrink!, false),
    (KenCarp4(), resize_shrink!, false),
    (TRBDF2(autodiff = AutoForwardDiff(chunksize = 2)), resize_shrink!, false),
    # With `alias_u0` the shrink also shrinks `prob.u0`.
    (Rodas5P(), resize_shrink!, true),
    (Rodas5P(), deleteat_shrink!, true),
    (Rodas5P(autodiff = AutoForwardDiff(chunksize = 5)), resize_shrink!, true),
    (TRBDF2(), resize_shrink!, true),
)

@testset "$(nameof(typeof(alg))) $spec $shrink! alias_u0=$alias_u0" for (alg, shrink!, alias_u0) in cases,
        spec in (SciMLBase.FullSpecialize, SciMLBase.AutoSpecialize)

    prob = ODEProblem{true, spec}(f!, collect(1.0:10.0), (0.0, 1.0))
    integ = init(
        prob, alg; abstol = 1.0e-10, reltol = 1.0e-10, dtmax = 1.0e-4, tstops = [0.5],
        alias = ODEAliasSpecifier(; alias_u0)
    )
    while integ.t < 0.5
        step!(integ)
    end
    shrink!(integ)
    SciMLBase.derivative_discontinuity!(integ, true)
    u_mid = copy(integ.u)
    solve!(integ)
    @test SciMLBase.successful_retcode(integ.sol)
    @test integ.t == 1.0

    ref = solve(
        ODEProblem{true, spec}(f!, u_mid, (0.5, 1.0)), alg;
        abstol = 1.0e-10, reltol = 1.0e-10, dtmax = 1.0e-4
    )
    @test integ.u ≈ ref.u[end] rtol = 1.0e-5
    exact = exp([-1.0 0.25; 0.25 -2.0] * 0.5) * u_mid
    @test integ.u ≈ exact rtol = 1.0e-7
    @test ref.u[end] ≈ exact rtol = 1.0e-7
end

@testset "Public Jacobian config pair" for backend in (
        AutoForwardDiff(), AutoForwardDiff(chunksize = 5), AutoFiniteDiff(),
    )
    integ = init(ODEProblem(f!, collect(1.0:10.0), (0.0, 1.0)), Rodas5P(; autodiff = backend))
    uf = SciMLBase.UJacobianWrapper(integ.f, integ.t, integ.p)
    du = zero(integ.u)
    config = OrdinaryDiffEqDifferentiation.build_jac_config(
        integ.alg, integ.f, uf, du, integ.uprev, integ.u, zero(du), zero(du)
    )
    @test config isa Tuple{Any, Any}
    @test length(config) == 2
    @test config[1] === config[2]
    J = zeros(10, 10)
    for prep in config[1:2]
        DI.jacobian!(uf, du, J, prep, integ.alg.autodiff, integ.u)
        @test diag(J) ≈ -collect(1.0:10.0)
    end
end
