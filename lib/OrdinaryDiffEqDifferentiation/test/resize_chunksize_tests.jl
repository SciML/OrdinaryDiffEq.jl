using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, SciMLBase, ADTypes, Test

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
        prob, alg; abstol = 1.0e-10, reltol = 1.0e-10, tstops = [0.5],
        alias = ODEAliasSpecifier(; alias_u0)
    )
    while integ.t < 0.5
        step!(integ)
    end
    shrink!(integ)
    u_mid = copy(integ.u)
    solve!(integ)
    @test SciMLBase.successful_retcode(integ.sol)
    @test integ.t == 1.0

    ref = solve(
        ODEProblem{true, spec}(f!, u_mid, (0.5, 1.0)), alg;
        abstol = 1.0e-10, reltol = 1.0e-10
    )
    # TRBDF2 itself is only accurate to about 3e-6 here, on a fresh solve as well.
    @test integ.u ≈ ref.u[end] rtol = 1.0e-5
end
