using OrdinaryDiffEqSDIRK
using Test

const SDIRK = OrdinaryDiffEqSDIRK

# A valid embedded pair has sum(bi) = sum(b̂i) = 1, hence sum(btilde) = 0.
adaptive_algs = [
    KenCarp3(), KenCarp4(), KenCarp5(), KenCarp47(), KenCarp58(),
    Kvaerno3(), Kvaerno4(), Kvaerno5(),
    ESDIRK54I8L2SA(), ESDIRK436L2SA2(), ESDIRK437L2SA(), ESDIRK547L2SA2(),
    ESDIRK659L2SA(),
    Hairer4(), Hairer42(), SDIRK2(), Cash4(), TRBDF2(),
]

@testset "Tableau embedded-pair consistency" begin
    for alg in adaptive_algs
        tab = SDIRK.ESDIRKIMEXTableau(alg, Float64, Float64)
        @testset "$(nameof(typeof(alg)))" begin
            @test isapprox(sum(tab.bi), 1; atol = 1e-9)
            if alg isa ESDIRK659L2SA
                @test_broken isapprox(sum(tab.btilde), 0; atol = 1e-9)  # #3659
            else
                @test isapprox(sum(tab.btilde), 0; atol = 1e-9)
            end
        end
    end
end
