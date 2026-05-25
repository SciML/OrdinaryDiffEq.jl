module EnzymeMakeSolutionZeroTests

using Test
using DiffEqBase
import ChainRulesCore, Enzyme  # triggers DiffEqBaseEnzymeExt
using SciMLBase

const EXT = Base.get_extension(DiffEqBase, :DiffEqBaseEnzymeExt)

@testset "EnzymeExt._make_solution_zero preserves prob.p / prob.u0 aliasing (NS#937)" begin
    # Regression for SciML/NonlinearSolve.jl#937 (ported to DiffEqBase). The
    # reverse rule builds the return-value shadow via `make_zero(sol)`.
    # A plain `Enzyme.make_zero` recursively zeros every mutable field of
    # `sol`, including `sol.prob.p` and `sol.prob.u0`, which the outer caller
    # has already registered as active shadows for the `p` / `u0` arguments.
    # Severing that aliasing means any cotangent written into the returned
    # `sol.prob.p` (or `.u0`) by a downstream consumer lands in a dangling
    # buffer instead of the one the outer Enzyme tape is tracking, silently
    # dropping that gradient contribution.
    #
    # `_make_solution_zero` pre-seeds the `make_zero` seen-set with identity
    # entries for `prob.p` and `prob.u0` so the recursion short-circuits and
    # the original buffers are reused verbatim in the shadow.

    f(du, u, p, t) = (du .= p .* u; nothing)
    u0 = [1.0, 1.0]
    p = [0.5, 0.25]
    tspan = (0.0, 1.0)
    prob = ODEProblem(f, u0, tspan, p)
    sol = SciMLBase.build_solution(
        prob, nothing, [0.0, 1.0], [copy(u0), copy(u0)];
        calculate_error = false,
    )

    # Naive `Enzyme.make_zero` allocates fresh buffers for prob.p / prob.u0.
    dsol_naive = Enzyme.make_zero(sol)
    @test objectid(dsol_naive.prob.p) != objectid(sol.prob.p)
    @test objectid(dsol_naive.prob.u0) != objectid(sol.prob.u0)

    # The extension helper keeps them aliased to the primal.
    dsol = EXT._make_solution_zero(sol)
    @test objectid(dsol.prob.p) == objectid(sol.prob.p)
    @test objectid(dsol.prob.u0) == objectid(sol.prob.u0)
    @test dsol.prob.p === sol.prob.p
    @test dsol.prob.u0 === sol.prob.u0
    # The actual derivative-carrying field (u) is still a fresh zero buffer.
    @test objectid(dsol.u) != objectid(sol.u)
    @test all(uu -> all(iszero, uu), dsol.u)

    # Guard the `nothing`/non-mutable path: a problem with `nothing` p
    # must not crash the pre-seed helper.
    prob_nop = ODEProblem(f, u0, tspan, SciMLBase.NullParameters())
    sol_nop = SciMLBase.build_solution(
        prob_nop, nothing, [0.0, 1.0], [copy(u0), copy(u0)];
        calculate_error = false,
    )
    dsol_nop = EXT._make_solution_zero(sol_nop)
    @test dsol_nop.prob.u0 === sol_nop.prob.u0
end

end
