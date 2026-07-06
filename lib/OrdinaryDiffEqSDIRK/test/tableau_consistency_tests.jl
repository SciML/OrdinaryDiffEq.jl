using OrdinaryDiffEqSDIRK
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
