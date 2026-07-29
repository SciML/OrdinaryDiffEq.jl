using OrdinaryDiffEqNonlinearSolve: default_nlsolve
using SciMLBase
using NonlinearSolve
using Test

# `default_nlsolve` selects the nonlinear algorithm the DAE-initialization OverrideInit path
# solves the init problem with. A `HomotopyProblem` (and an `SCCNonlinearProblem` that
# contains `HomotopyProblem` blocks) must be initialized by continuation.

Hf(u, p, λ) = u .^ 2 .- 4λ
hp = HomotopyProblem(Hf, [1.0]; λspan = (0.0, 1.0))
np = NonlinearProblem((u, p) -> u .^ 2 .- 4, [1.5])

@testset "HomotopyProblem -> continuation" begin
    for iip in (Val(true), Val(false))
        @test default_nlsolve(nothing, iip, [1.0], hp) isa HomotopyPolyAlgorithm
    end
end

@testset "inner corrector AD is threaded (ForwardDiff-hostile residual)" begin
    # a residual that errors on dual numbers can only be finite-differenced
    badf(x::Float64) = atan(x - 3)
    badf(x) = throw(ArgumentError("residual seen a dual number (ForwardDiff)"))
    H(u, p, λ) = [(1 - λ) * u[1] + λ * badf(u[1])]
    hp2 = HomotopyProblem(H, [12.0]; λspan = (0.0, 1.0))

    sol_fd = solve(hp2, default_nlsolve(nothing, Val(false), [12.0], hp2, false); abstol = 1.0e-10)
    @test SciMLBase.successful_retcode(sol_fd)
    @test sol_fd.u[1] ≈ 3.0 atol = 1.0e-6

    threw = false
    try
        solve(hp2, default_nlsolve(nothing, Val(false), [12.0], hp2, true); abstol = 1.0e-10)
    catch
        threw = true
    end
    @test threw
end

@testset "plain NonlinearProblem keeps its Newton polyalgorithm default" begin
    @test default_nlsolve(nothing, Val(false), [1.5], np) !== nothing
    @test default_nlsolve(nothing, Val(true), [1.5], np) !== nothing
    @test !(default_nlsolve(nothing, Val(false), [1.5], np) isa HomotopyPolyAlgorithm)
end

@testset "SCCNonlinearProblem: nothing iff it has a HomotopyProblem block" begin
    explicit = (Returns(nothing), Returns(nothing))
    scc_hom = SCCNonlinearProblem((np, hp), explicit, nothing, Val(true))
    scc_plain = SCCNonlinearProblem((np, np), explicit, nothing, Val(true))
    u = [1.5, 1.0]
    @test default_nlsolve(nothing, Val(false), u, scc_hom) === nothing
    @test default_nlsolve(nothing, Val(true), u, scc_hom) === nothing
    @test default_nlsolve(nothing, Val(false), u, scc_plain) !== nothing
    @test default_nlsolve(nothing, Val(true), u, scc_plain) !== nothing
end

@testset "stateless init problem selects no solver" begin
    # `u === nothing` (an `initializeprob` with `state_values(...) === nothing`): the
    # OrdinaryDiffEqCore method covering it is more specific in `u` while these are more
    # specific in the problem type, so without an explicit intersection method the call is
    # an ambiguity error rather than the `nothing` Core answers everywhere else.
    scc = SCCNonlinearProblem((np, hp), (Returns(nothing), Returns(nothing)), nothing, Val(true))
    for iip in (Val(true), Val(false))
        @test default_nlsolve(nothing, iip, nothing, hp) === nothing
        @test default_nlsolve(nothing, iip, nothing, scc) === nothing
    end
end
