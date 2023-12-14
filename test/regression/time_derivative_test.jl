using OrdinaryDiffEq, StaticArrays, Test

function time_derivative(du, u, p, t)
    du[1] = -t
end
function time_derivative_static(u, p, t)
    SVector(-t)
end
function time_derivative_analytic(u0, p, t)
    u0 .- t .^ 2 ./ 2
end

const CACHE_TEST_ALGS = [Euler(), Midpoint(), RK4(), SSPRK22(), SSPRK33(), SSPRK53(),
    SSPRK63(), SSPRK73(), SSPRK83(), SSPRK43(), SSPRK432(), SSPRK932(), SSPRK54(),
    SSPRK104(), CarpenterKennedy2N54(),
    BS3(), BS5(), DP5(), DP8(), Feagin10(), Feagin12(), Feagin14(), TanYam7(),
    Tsit5(), TsitPap8(), Vern6(), Vern7(), Vern8(), Vern9(), OwrenZen3(), OwrenZen4(),
    OwrenZen5()]

tspan = (0.0, 1.0)

for (ff_time_derivative, u0) in ((ODEFunction(time_derivative,
            analytic = time_derivative_analytic), [1.0]),
    (ODEFunction(time_derivative_static,
            analytic = time_derivative_analytic),
        SVector(1.0)))
    @info "StaticArrays?: $(u0 isa StaticArray)"

    prob = ODEProblem(ff_time_derivative, u0, tspan)

    for _autodiff in (true, false)
        @info "autodiff=$(_autodiff)"

        @show Rosenbrock23
        sol = solve(prob, Rosenbrock32(autodiff = _autodiff), reltol = 1e-9, abstol = 1e-9)
        @test sol.errors[:final] < 1e-5 * 10^(!_autodiff)
        sol = solve(prob, Rosenbrock23(autodiff = _autodiff), reltol = 1e-9, abstol = 1e-9)
        @test sol.errors[:final] < 1e-10 * 10_000^(!_autodiff)

        @show Rodas4
        sol = solve(prob, Rodas4(autodiff = _autodiff), reltol = 1e-9, abstol = 1e-9)
        @test sol.errors[:final] < 1e-10 * 100_000^(!_autodiff)

        @show Rodas5
        sol = solve(prob, Rodas5(autodiff = _autodiff), reltol = 1e-9, abstol = 1e-9)
        @test sol.errors[:final] < 1e-10 * 1_000_000_000^(!_autodiff)

        @show Veldd4
        sol = solve(prob, Veldd4(autodiff = _autodiff), reltol = 1e-9, abstol = 1e-9)
        @test sol.errors[:final] < 1e-10 * 100_000^(!_autodiff)

        @show KenCarp3
        sol = solve(prob, KenCarp3(autodiff = _autodiff), reltol = 1e-12, abstol = 1e-12)
        @test length(sol) > 2
        @test sol.errors[:final] < 1e-10

        @show KenCarp4
        sol = solve(prob, KenCarp4(autodiff = _autodiff), reltol = 1e-12, abstol = 1e-12)
        @test length(sol) > 2
        @test sol.errors[:final] < 1e-10

        @show KenCarp47
        sol = solve(prob, KenCarp47(autodiff = _autodiff), reltol = 1e-12, abstol = 1e-12)
        @test length(sol) > 2
        @test sol.errors[:final] < 1e-10

        @show KenCarp5
        sol = solve(prob, KenCarp5(autodiff = _autodiff), reltol = 1e-12, abstol = 1e-12)
        @test length(sol) > 2
        @test sol.errors[:final] < 1e-10

        @show KenCarp58
        sol = solve(prob, KenCarp58(autodiff = _autodiff), reltol = 1e-12, abstol = 1e-12)
        @test length(sol) > 2
        @test sol.errors[:final] < 1e-10

        @show TRBDF2
        sol = solve(prob, TRBDF2(autodiff = _autodiff), reltol = 1e-9, abstol = 1e-9)
        @test sol.errors[:final] < 1e-10

        @show ImplicitEuler
        sol = solve(prob, ImplicitEuler(autodiff = _autodiff), dt = 1 / 10)
        @test sol.errors[:final] < 1e-1

        @show Trapezoid
        sol = solve(prob, Trapezoid(autodiff = _autodiff), dt = 1 / 10)
        @test sol.errors[:final] < 1e-12
    end

    @show Euler
    sol = solve(prob, Euler(), dt = 1 / 100)
    @test sol.errors[:final] < 6e-3

    for alg in CACHE_TEST_ALGS
        @show alg
        sol = solve(prob, alg, dt = 1 / 10)
        if !(alg isa Euler)
            @test sol.errors[:final] < 4e-14
        end
    end
end
