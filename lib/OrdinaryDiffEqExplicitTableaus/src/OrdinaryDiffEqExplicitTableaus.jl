module OrdinaryDiffEqExplicitTableaus

import DiffEqBase

export Baker10, BogakiShampine3, BogakiShampine5, Butcher6, Butcher62, Butcher63,
    Butcher7, CashKarp, Cassity5, Chummund6, Chummund62, ClassicVerner6, ClassicVerner7,
    ClassicVerner8, CooperVerner8, CooperVerner82, Curtis10, Curtis8,
    DormandLockyerMcCorriganPrince6, DormandPrince6, DormandPrince8, DormandPrince8_64bit,
    Dverk, dverk78, EnrightVerner7, EnrightVerner8, Euler, Feagin10, Feagin12, Feagin14,
    Hairer10, Heun, Huta6, Huta62, Kutta3, Lawson5, Lawson6, LobattoIIICStar2,
    LutherKonen5, LutherKonen52, LutherKonen53, MikkawyEisa, Ono10, Ono12, Papakostas6,
    PapakostasPapaGeorgiou5, PapakostasPapaGeorgiou52, Ralston, Ralston4, RK4, RK438Rule,
    RKF4, RKF5, RKF8, RKO65, RungeFirst5, Sharp9, SharpSmart5, SharpSmart7, SharpVerner6,
    SharpVerner7, SSPRK104, SSPRK22, SSPRK33, SSPRK43,
    TanakaKasugaYamashitaYazaki6A, TanakaKasugaYamashitaYazaki6B,
    TanakaKasugaYamashitaYazaki6C, TanakaKasugaYamashitaYazaki6D,
    TanakaYamashitaEfficient7, TanakaYamashitaStable7, Tsitouras5, Tsitouras9, Tsitouras92,
    TsitourasPapakostas6, TsitourasPapakostas8, Verner6, Verner7, Verner8, Verner916,
    Verner9162, VernerEfficient6, VernerEfficient7, VernerEfficient9, VernerRobust6,
    VernerRobust7, VernerRobust9

include("tableaus_low_order.jl")   # Orders 1-4: Euler, Heun, RK4, SSPRK, etc.
include("tableaus_order5.jl")      # Order 5: RKF5, Tsitouras5, BogakiShampine5, CashKarp, etc.
include("tableaus_order6.jl")      # Order 6: Butcher6, Verner6, DormandPrince6, etc.
include("tableaus_order7.jl")      # Order 7: Verner7, EnrightVerner7, SharpSmart7, etc.
include("tableaus_order8_9.jl")    # Orders 8-9: CooperVerner8, DormandPrince8, Verner9, etc.
include("tableaus_high_order.jl")  # Orders 10-14: Feagin10/12/14, Baker10, Hairer10, etc.
include("tableaus_classic.jl")     # Classic wrappers: Butcher7, Dverk, ClassicVerner6/7/8, RKF8, etc.

end
