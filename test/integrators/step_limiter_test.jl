using OrdinaryDiffEq, Test

# define the counting variable
const STEP_LIMITER_VAR = Ref(0)
# define the step_limiter! function which just counts the number of step_limiter calls
step_limiter!(u, integrator, p, t) = STEP_LIMITER_VAR[] += 1

# This function tests the step limiter functionality of an ODE solver.
function test_step_limiter(alg_type)
    STEP_LIMITER_VAR[] = 0 # reset the counting variable
    prob = ODEProblem((du, u, p, t) -> du .= u, [1.0], (0.0, 1.0))

    sol = solve(prob, alg_type(; step_limiter!), dt = 0.1)

    n = sol.stats.naccept + sol.stats.nreject
    @test n == STEP_LIMITER_VAR[]
end

@testset "Step_limiter Test" begin
    # it only catches the most basic errors, i.e. if the step_limiter! function is not called
    # or called more then one time

    # test the step_limiter! function 
    alg_types = [
        QNDF1, QNDF2, QNDF, FBDF, ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2,
        SDIRK2, SDIRK22, ABDF2, Feagin10, Feagin12, Feagin14,
        KenCarp3, KenCarp4, KenCarp5, Kvaerno3, Kvaerno4, Kvaerno5,
        Rosenbrock23, Rosenbrock32, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42,
        Rodas4P, Rodas4P2, Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr, 
        AdaptiveRadau, RadauIIA9, RadauIIA5, RadauIIA3, SIR54, Alshina2, Alshina3, Heun, Ralston, Midpoint, RK4,
        OwrenZen3, OwrenZen4, OwrenZen5,
        BS3, DP5, Tsit5, DP8, TanYam7, TsitPap8, FRK65, PFRK87, BS5, Vern6, Vern7,
        Vern8, Vern9, QPRK98, SSPRKMSVS43, SSPRKMSVS32, SSPRK432, SSPRK43,
        RDPK3SpFSAL35, RDPK3Sp35, NDBLSRK124, NDBLSRK134, DGLDDRK73_C,
        DGLDDRK84_C, DGLDDRK84_F, SHLDDRK64, RDPK3Sp49, RDPK3SpFSAL49, RDPK3Sp510, RDPK3SpFSAL510,
        Alshina6, RKM, MSRK5, MSRK6, Anas5, RKO65, RK46NL, ORK256, KYK2014DGSSPRK_3S2,
        SSPRK22, SSPRK104, SSPRK54, SSPRK932, SSPRK83, SSPRK73, SSPRK63, SSPRK53_H,
        SSPRK53_2N2, SSPRK53_2N1, SSPRK53, SSPRK33, SHLDDRK_2N, SHLDDRK52, KYKSSPRK42,
        CarpenterKennedy2N54, CFRLDDRK64, TSLDDRK74, ParsaniKetchesonDeconinck3S32,
        ParsaniKetchesonDeconinck3S82,
        ParsaniKetchesonDeconinck3S53, ParsaniKetchesonDeconinck3S173, ParsaniKetchesonDeconinck3S94,
        ParsaniKetchesonDeconinck3S184, ParsaniKetchesonDeconinck3S105, ParsaniKetchesonDeconinck3S205,
        CKLLSRK43_2, CKLLSRK54_3C, CKLLSRK95_4S, CKLLSRK95_4C, CKLLSRK95_4M, CKLLSRK54_3C_3R,
        CKLLSRK54_3M_3R, CKLLSRK54_3N_3R, CKLLSRK85_4C_3R, CKLLSRK85_4M_3R, CKLLSRK85_4P_3R,
        CKLLSRK54_3N_4R, CKLLSRK54_3M_4R, CKLLSRK65_4M_4R, CKLLSRK85_4FM_4R, CKLLSRK75_4M_5R] #Stepanov5    

    for alg_type in alg_types
        test_step_limiter(alg_type)
    end
end
