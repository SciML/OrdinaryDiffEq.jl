module OrdinaryDiffEqExplicitTableaus

import DiffEqBase

include("tableaus_low_order.jl")   # Orders 1-4: Euler, Heun, RK4, SSPRK, etc.
include("tableaus_order5.jl")      # Order 5: RKF5, Tsitouras5, BogakiShampine5, CashKarp, etc.
include("tableaus_order6.jl")      # Order 6: Butcher6, Verner6, DormandPrince6, etc.
include("tableaus_order7.jl")      # Order 7: Verner7, EnrightVerner7, SharpSmart7, etc.
include("tableaus_order8_9.jl")    # Orders 8-9: CooperVerner8, DormandPrince8, Verner9, etc.
include("tableaus_high_order.jl")  # Orders 10-14: Feagin10/12/14, Baker10, Hairer10, etc.
include("tableaus_classic.jl")     # Classic wrappers: Butcher7, Dverk, ClassicVerner6/7/8, RKF8, etc.

end
