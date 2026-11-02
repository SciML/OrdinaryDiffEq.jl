using OrdinaryDiffEqSDIRK, LinearAlgebra
using Test

const SDIRK = OrdinaryDiffEqSDIRK

all_algs = [
    ImplicitEuler(), ImplicitMidpoint(), Trapezoid(), TRBDF2(), SDIRK2(),
    Kvaerno3(), KenCarp3(), Cash4(), Hairer4(), Hairer42(), SSPSDIRK2(),
    Kvaerno4(), Kvaerno5(), KenCarp4(), KenCarp47(), KenCarp5(), KenCarp58(),
    ESDIRK54I8L2SA(), SFSDIRK4(), SFSDIRK5(), SFSDIRK6(), SFSDIRK7(), SFSDIRK8(),
    ESDIRK325L2SA(), ESDIRK436L2SA2(), ESDIRK437L2SA(), ESDIRK547L2SA2(), ESDIRK659L2SA(),
    CFNLIRK3(), ARS343(), ARS222(), ARS232(), ARS443(),
    IMEXSSP222(), IMEXSSP2322(), IMEXSSP3332(), IMEXSSP3433(), BHR553(),
]

@testset "Tableau embedded-pair consistency" begin
    for alg in all_algs
        tab = SDIRK.ESDIRKIMEXTableau(alg, Float64, Float64)
        @testset "$(nameof(typeof(alg)))" begin
            @test isapprox(sum(tab.bi), 1; atol = 1.0e-9)
            isempty(tab.btilde) && continue
            @test isapprox(sum(tab.btilde), 0; atol = 1.0e-9)
        end
    end
end

@testset "A step ends with the derivative at the new solution" begin
    forced!(du, u, p, t) = (du[1] = cos(t) - u[1]^2; nothing)
    forced(u, p, t) = [cos(t) - u[1]^2]
    cases = Any[(alg, forced!) for alg in all_algs]
    append!(cases, [(alg, forced) for alg in (SDIRK2(), SSPSDIRK2(), SFSDIRK4())])
    for (alg, rhs) in cases
        @testset "$(nameof(typeof(alg)))" begin
            integrator = init(ODEProblem(rhs, [1.0], (0.0, 1.0)), alg; dt = 0.1, adaptive = false)
            step!(integrator)
            @test integrator.fsallast ≈ forced(integrator.u, nothing, integrator.t)
        end
    end
    # `2u' = -2u`: storing `f(u)` instead of `du/dt` drops the dense output to first order.
    decay!(du, u, p, t) = (du[1] = -2 * u[1]; nothing)
    decay(u, p, t) = [-2 * u[1]]
    mm_algs = (
        SFSDIRK4(), SFSDIRK5(), SFSDIRK6(), SFSDIRK7(), SFSDIRK8(),
        SDIRK2(), SSPSDIRK2(), ImplicitMidpoint(),
    )
    for alg in mm_algs, rhs in (decay!, decay)
        @testset "$(nameof(typeof(alg))) mass matrix, $(nameof(rhs))" begin
            errs = map((1 / 32, 1 / 64)) do dt
                f = ODEFunction(rhs; mass_matrix = Diagonal([2.0]))
                sol = solve(ODEProblem(f, [1.0], (0.0, 1.0)), alg; dt, adaptive = false)
                maximum(t -> abs(sol(t)[1] - exp(-t)), (0.5:dt:(1 - dt)) .+ dt / 4)
            end
            @test log2(errs[1] / errs[2]) > 1.5
        end
    end
end

@testset "Non-IMEX methods solve a SplitODEProblem as f1 + f2" begin
    # A method with no explicit part has nowhere to put `f2`, so it integrates the combined
    # right-hand side, the way a non-split solve of the same method does. Before #4149 the
    # explicit accumulator ran anyway and indexed the empty `Ae`/`be`/`split_guess`.
    f1 = (du, u, p, t) -> (du .= -u)
    f2 = (du, u, p, t) -> (du .= 0.1 .* u)
    psplit = SplitODEProblem(f1, f2, [1.0], (0.0, 1.0))
    pcomb = ODEProblem((du, u, p, t) -> (du .= -0.9 .* u), [1.0], (0.0, 1.0))

    for alg in all_algs
        SDIRK.issplit(alg) && continue
        @testset "$(nameof(typeof(alg)))" begin
            a = solve(psplit, alg; dt = 0.01, adaptive = false)
            b = solve(pcomb, alg; dt = 0.01, adaptive = false)
            @test SciMLBase.successful_retcode(a)
            @test a.u[end] ≈ b.u[end] rtol = 1.0e-12
        end
    end
end
