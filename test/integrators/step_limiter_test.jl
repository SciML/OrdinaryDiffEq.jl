using OrdinaryDiffEq, Test

# define the counting variable
const STEP_LIMITER_VAR = Ref(0)
# define the step_limiter! function which just counts the number of step_limiter calls
step_limiter!(u, integrator, p, t) = STEP_LIMITER_VAR[] += 1


# This function tests the step limiter functionality of an ODE solver.
function test_step_limiter(alg_type, adaptive = true)
    STEP_LIMITER_VAR[] = 0 # reset the counting variable
    prob = ODEProblem((du, u, p, t) -> du .= u, [1.0], (0.0, 5.0))

    if adaptive == true
        sol = solve(prob, alg_type(; step_limiter!))
    else
        sol = solve(prob, alg_type(; step_limiter!), dt = 1)
    end

    n = sol.stats.naccept + sol.stats.nreject
    @test n == STEP_LIMITER_VAR[]
end


@testset "Step_limiter Test" begin
    # it only catches the most basic errors, i.e. if the step_limiter! function is not called
    # or called more then one time

    # I think all the algorithms with step_limiter! function are tested here
    # (besides Low Storage Runge-Kutta methods)

    STEP_LIMITER_VAR[] = 0
    
    # test the step_limiter! function for adaptive algorithms
    adaptive = [RadauIIA5, RadauIIA3, SIR54, Alshina2, Alshina3, Heun, Ralston, Midpoint, RK4,
                OwrenZen3, OwrenZen4, OwrenZen5, 
                BS3, DP5, Tsit5, DP8, TanYam7, TsitPap8, FRK65, PFRK87, BS5, Vern6, Vern7,
                Vern8, Vern9, QPRK98, SSPRKMSVS43, SSPRKMSVS32, SSPRK432, SSPRK43]
    for alg_type in adaptive
        test_step_limiter(alg_type)
    end


    # test the step_limiter! function for non adaptive algorithms
    non_adaptive = [Alshina6, RKM, MSRK5, MSRK6, Anas5, RKO65, RK46NL, ORK256, KYK2014DGSSPRK_3S2,
                    SSPRK22, SSPRK104, SSPRK54, SSPRK932, SSPRK83, SSPRK73, SSPRK63, SSPRK53_H, SSPRK53_2N2,
                    SSPRK53_2N1, SSPRK53, SSPRK33, SHLDDRK_2N, SHLDDRK52, KYKSSPRK42]
    for alg_type in non_adaptive
        test_step_limiter(alg_type, false)
    end

    #=  
        maybe it would be nicer to check for adaptive and non adaptive in the 
        test_step_limiter function, but I don't know how to do that / couldn't make 
        "isadaptive()" work 
    =#
end
