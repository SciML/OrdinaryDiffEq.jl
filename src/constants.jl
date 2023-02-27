"""
constructDormandPrince()

Constructs the tableau object for the Dormand-Prince Order 4/5 method.
"""
function constructDormandPrince(T::Type = Float64)
    A = [0 0 0 0 0 0 0
         1//5 0 0 0 0 0 0
         3//40 9//40 0 0 0 0 0
         44//45 -56//15 32//9 0 0 0 0
         19372//6561 -25360//2187 64448//6561 -212//729 0 0 0
         9017//3168 -355//33 46732//5247 49//176 -5103//18656 0 0
         35//384 0 500//1113 125//192 -2187//6784 11//84 0]
    c = [0; 1 // 5; 3 // 10; 4 // 5; 8 // 9; 1; 1]
    α = [35 // 384; 0; 500 // 1113; 125 // 192; -2187 // 6784; 11 // 84; 0]
    αEEst = [5179 // 57600; 0; 7571 // 16695; 393 // 640; -92097 // 339200; 187 // 2100;
             1 // 40]
    A = map(T, A)
    α = map(T, α)
    αEEst = map(T, αEEst)
    c = map(T, c)
    return (DiffEqBase.ExplicitRKTableau(A, c, α, 5, αEEst = αEEst, adaptiveorder = 4,
                                         fsal = true, stability_size = 3.3066))
end

"""
ODE_DEFAULT_TABLEAU

Sets the default tableau for the ODE solver. Currently Dormand-Prince 4/5.
"""
const ODE_DEFAULT_TABLEAU = constructDormandPrince()
