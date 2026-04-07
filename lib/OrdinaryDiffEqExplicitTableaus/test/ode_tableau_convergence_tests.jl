using OrdinaryDiffEq, DiffEqDevTools, Test, Random
using ODEProblemLibrary: prob_ode_2Dlinear, prob_ode_linear,
    prob_ode_bigfloatlinear, prob_ode_bigfloat2Dlinear

const ET = OrdinaryDiffEqExplicitTableaus

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

for tab in [ET.Ralston4(Rational{BigInt}), ET.TsitourasPapakostas6(Rational{BigInt}), ET.DormandLockyerMcCorriganPrince6(Rational{BigInt})]
    @test all(i -> residual_order_condition(tab, i, +, abs) < 10eps(1.0), 1:(tab.order))
    if tab.adaptiveorder != 0
        @test all(
            i -> residual_order_condition(tab, i, +, abs; embedded = true) < 10eps(1.0),
            tab.adaptiveorder)
    end
end

for i in 1:2 # 1 = num, 2 = ExplicitRK
    global dts
    if i > 1
        prob = probArr[2]
        bigprob = bigprobArr[2]
    else
        prob = probArr[1]
        bigprob = bigprobArr[1]
    end

    # Order 2

    tabalg = ExplicitRK(tableau = ET.Heun())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 2) < testTol

    tabalg = ExplicitRK(tableau = ET.Ralston())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 2) < testTol

    tabalg = ExplicitRK(tableau = ET.SSPRK22())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 2) < testTol

    # Order 3

    tabalg = ExplicitRK(tableau = ET.BogakiShampine3())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 3) < testTol

    tabalg = ExplicitRK(tableau = ET.SSPRK33())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 3) < testTol

    tabalg = ExplicitRK(tableau = ET.SSPRK43())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 3) < testTol

    # Order 4

    tabalg = ExplicitRK(tableau = ET.RKF4())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 4) < testTol

    tabalg = ExplicitRK(tableau = ET.Ralston4())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 4) < testTol

    tabalg = ExplicitRK(tableau = ET.SSPRK104())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 4) < testTol

    # Order 5

    dts = 1 .// 2 .^ (7:-1:4)
    tabalg = ExplicitRK(tableau = ET.RKF5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    dts = 1 .// 2 .^ (7:-1:4)
    tabalg = ExplicitRK(tableau = ET.DormandPrince())
    sim = test_convergence(dts, prob, tabalg)
    sim2 = test_convergence(dts, prob, DP5())
    @test (abs(sim.𝒪est[:l∞] - 5) < testTol && (maximum(sim[end][end] - sim2[end][end]) < 1e-10))

    tabalg = ExplicitRK(tableau = ET.CashKarp())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    dts = 1 .// 2 .^ (7:-1:4)
    tabalg = ExplicitRK(tableau = ET.RungeFirst5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.Cassity5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.Lawson5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.LutherKonen5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.LutherKonen52())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.LutherKonen53())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.PapakostasPapaGeorgiou5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.PapakostasPapaGeorgiou52())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    dts = 1 .// 2 .^ (6:-1:4)
    tabalg = ExplicitRK(tableau = ET.Tsitouras5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    dts = 1 .// 2 .^ (6:-1:4)
    tabalg = ExplicitRK(tableau = ET.BogakiShampine5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.SharpSmart5())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    tabalg = ExplicitRK(tableau = ET.RKO65())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5) < testTol

    # Order 6

    dts = 1 .// 2 .^ (6:-1:4)
    tabalg = ExplicitRK(tableau = ET.Butcher6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.Butcher62())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    dts = 1 .// 2 .^ (6:-1:4)
    tabalg = ExplicitRK(tableau = ET.Butcher63())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    dts = 1 .// 2 .^ (5:-1:1)
    tabalg = ExplicitRK(tableau = ET.DormandPrince6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7) < testTol + 0.1 # Better on linear

    tabalg = ExplicitRK(tableau = ET.SharpVerner6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.Verner916())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.Verner9162())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.VernerRobust6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.VernerEfficient6(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6.6) < testTol

    tabalg = ExplicitRK(tableau = ET.Papakostas6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.Lawson6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    dts = 1 .// 2 .^ (3:-1:1)
    tabalg = ExplicitRK(tableau = ET.TsitourasPapakostas6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6.7) < testTol # Better on linear

    dts = 1 .// 2 .^ (5:-1:1)
    tabalg = ExplicitRK(tableau = ET.DormandLockyerMcCorriganPrince6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.TanakaKasugaYamashitaYazaki6D())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.TanakaKasugaYamashitaYazaki6C())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.TanakaKasugaYamashitaYazaki6B())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.TanakaKasugaYamashitaYazaki6A())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.MikkawyEisa())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6.53) < testTol # Odd behavior

    tabalg = ExplicitRK(tableau = ET.Chummund6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.Chummund62())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.Huta6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 5.5) < testTol # Low convergence, error noted in Stone notes

    dts = 1 .// 2 .^ (5:-1:1)
    tabalg = ExplicitRK(tableau = ET.Huta62())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.Verner6())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6.7) < testTol # Better on linear

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.Dverk())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    tabalg = ExplicitRK(tableau = ET.ClassicVerner6())
    sim = test_convergence(dts, prob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6) < testTol

    # Order 7

    dts = 1 .// 2 .^ (5:-1:1)
    tabalg = ExplicitRK(tableau = ET.Butcher7())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7) < testTol

    dts = 1 .// 2 .^ (5:-1:2)
    tabalg = ExplicitRK(tableau = ET.ClassicVerner7())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7) < testTol

    tabalg = ExplicitRK(tableau = ET.VernerRobust7())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7) < testTol

    dts = 1 .// 2 .^ (5:-1:1)
    tabalg = ExplicitRK(tableau = ET.EnrightVerner7())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7.15) < testTol # Better on linear

    tabalg = ExplicitRK(tableau = ET.TanakaYamashitaStable7())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7.3) < testTol

    tabalg = ExplicitRK(tableau = ET.TanakaYamashitaEfficient7(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7) < testTol

    dts = 1 .// 2 .^ (8:-1:3)
    tabalg = ExplicitRK(tableau = ET.SharpSmart7(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7) < testTol

    dts = 1 .// 2 .^ (3:-1:1)
    tabalg = ExplicitRK(tableau = ET.SharpVerner7())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 6.5) < testTol # Coefficients aren't accurate enough, drop off error

    tabalg = ExplicitRK(tableau = ET.VernerEfficient7())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7) < testTol

    tabalg = ExplicitRK(tableau = ET.Verner7())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 7) < testTol

    # Order 8
    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.ClassicVerner8())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    dts = 1 .// 2 .^ (5:-1:1)
    tabalg = ExplicitRK(tableau = ET.Verner8(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.CooperVerner8())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.CooperVerner82())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    tabalg = ExplicitRK(tableau = ET.TsitourasPapakostas8(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.dverk78(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.EnrightVerner8(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.Curtis8())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    dts = 1 .// 2 .^ (4:-1:1)
    tabalg = ExplicitRK(tableau = ET.RKF8(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8) < testTol

    tabalg = ExplicitRK(tableau = ET.DormandPrince8(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8.4) < testTol

    dts = 1 .// 2 .^ (3:-1:1)
    tabalg = ExplicitRK(tableau = ET.DormandPrince8_64bit(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 8.4) < testTol

    # Order 9

    tabalg = ExplicitRK(tableau = ET.VernerRobust9(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 9) < testTol

    tabalg = ExplicitRK(tableau = ET.VernerEfficient9(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 9) < testTol

    dts = 1 .// 2 .^ (3:-1:1)
    tabalg = ExplicitRK(tableau = ET.Sharp9(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 9) < testTol

    dts = 1 .// 2 .^ (2:-1:1)
    tabalg = ExplicitRK(tableau = ET.Tsitouras9(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 10.5) < testTol #Only works to Float64

    dts = 1 .// 2 .^ (3:-1:1)
    tabalg = ExplicitRK(tableau = ET.Tsitouras92(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 9) < testTol  #Only works to Float64

    ## Order 10

    dts = 1 .// 2 .^ (5:-1:1)
    tabalg = ExplicitRK(tableau = ET.Curtis10())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 10) < testTol

    tabalg = ExplicitRK(tableau = ET.Ono10())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 10) < testTol

    dts = 1 .// 2 .^ (5:-1:1)
    tabalg = ExplicitRK(tableau = ET.Feagin10(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 10) < testTol

    tabalg = ExplicitRK(tableau = ET.Curtis10())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 10) < testTol

    dts = 1 .// 2 .^ (6:-1:1)
    tabalg = ExplicitRK(tableau = ET.Baker10())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 10.8) < testTol

    tabalg = ExplicitRK(tableau = ET.Hairer10())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 10.8) < testTol

    ## Order 12

    dts = 1 .// 2 .^ (6:-1:1)
    tabalg = ExplicitRK(tableau = ET.Feagin12(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 12.6) < testTol

    tabalg = ExplicitRK(tableau = ET.Ono12())
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 11.6) < testTol

    ## Order 14

    dts = 1 .// 2 .^ (6:-1:1)
    tabalg = ExplicitRK(tableau = ET.Feagin14(BigFloat))
    sim = test_convergence(dts, bigprob, tabalg)
    @test abs(sim.𝒪est[:l∞] - 15.5) < testTol
end
