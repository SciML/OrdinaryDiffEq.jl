using OrdinaryDiffEq, DiffEqDevTools, Test, Random
using ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear,
                         prob_ode_bigfloatlinear,
                         prob_ode_bigfloat2Dlinear

probArr = Vector{ODEProblem}(undef, 2)
bigprobArr = Vector{ODEProblem}(undef, 2)

probArr[1] = prob_ode_linear
probArr[2] = prob_ode_2Dlinear
bigprobArr[1] = prob_ode_bigfloatlinear
bigprobArr[2] = prob_ode_bigfloat2Dlinear
setprecision(400)
Random.seed!(100)
## Convergence Testing
println("Convergence Test on Linear")
dts = 1 .// 2 .^ (8:-1:4)
testTol = 0.3
superduperbool = Vector{Bool}(undef, 2)

for constructfun in filter(x -> startswith(string(x), "construct"), names(DiffEqDevTools))
    tab = getproperty(DiffEqDevTools, constructfun)(BigFloat)
    if tab.order < 12
        if constructfun in (:constructTsitouras9,  # order 1 ???!!!
            :constructTsitouras92)
            @info "Failed $constructfun..."
            @test_broken check_tableau(tab)
        else
            @info "Testing $constructfun..."
            @test check_tableau(tab)
        end
    end
end

# numerically test RKs with order >= 12
for i in 1:2 # 1 = num, 2 = ExplicitRK
    global dts
    if i > 1
        prob = probArr[2]
        bigprob = bigprobArr[2]
    else
        prob = probArr[1]
        bigprob = bigprobArr[1]
    end

    for tab in [
        constructFeagin12(BigFloat),
        constructFeagin14(BigFloat),
        constructOno12(BigFloat)
    ]
        @info "Testing..."
        tabalg = ExplicitRK(tableau = tab)
        sim = test_convergence(dts, bigprob, tabalg)
        @test sim.𝒪est[:l∞] >= tab.order || abs(sim.𝒪est[:l∞] - tab.order) < testTol
    end
end
