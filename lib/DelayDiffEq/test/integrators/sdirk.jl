using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve: NLNewton, NLFunctional, NLAnderson
using Test

const prob_ip = prob_dde_constant_1delay_ip
const prob_scalar = prob_dde_constant_1delay_scalar
const ts = 0:0.1:10

# ODE algorithms
noreuse = NLNewton(fast_convergence_cutoff = 0)
const working_algs = [
    ImplicitMidpoint(), SSPSDIRK2(), KenCarp5(nlsolve = noreuse),
    ImplicitEuler(nlsolve = noreuse), Trapezoid(nlsolve = noreuse),
    TRBDF2(nlsolve = noreuse), SDIRK2(nlsolve = noreuse),
    Kvaerno3(nlsolve = noreuse), KenCarp3(nlsolve = noreuse),
    Cash4(nlsolve = noreuse), Hairer4(nlsolve = noreuse),
    Hairer42(nlsolve = noreuse), Kvaerno4(nlsolve = noreuse), KenCarp4(nlsolve = noreuse),
    Kvaerno5(nlsolve = noreuse),
]

const broken_algs = [
    ImplicitEuler(), Trapezoid(),
    TRBDF2(), SDIRK2(),
    Kvaerno3(), KenCarp3(),
    Cash4(), Hairer4(), Hairer42(), Kvaerno4(), KenCarp4(),
    Kvaerno5(),
]

@testset "Algorithm $(nameof(typeof(alg)))" for alg in working_algs
    println(nameof(typeof(alg)))

    stepsalg = MethodOfSteps(alg)
    sol_ip = solve(prob_ip, stepsalg; dt = 0.1)
    sol_scalar = solve(prob_scalar, stepsalg; dt = 0.1)

    @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts)
    @test sol_ip.t ≈ sol_scalar.t
    @test sol_ip[1, :] ≈ sol_scalar.u
end

@testset "Algorithm $(nameof(typeof(alg)))" for alg in broken_algs
    println(nameof(typeof(alg)))

    stepsalg = MethodOfSteps(alg)
    sol_ip = solve(prob_ip, stepsalg; dt = 0.1)
    sol_scalar = solve(prob_scalar, stepsalg; dt = 0.1)

    @test_broken sol_ip(ts, idxs = 1) ≈ sol_scalar(ts)
    # this test is not broken for KenCarp4
    #if alg isa KenCarp4
    #  @test sol_ip.t ≈ sol_scalar.t
    #else
    @test_broken sol_ip.t ≈ sol_scalar.t
    #end
    @test_broken sol_ip[1, :] ≈ sol_scalar.u
end
