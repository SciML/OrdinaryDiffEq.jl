using DiffEqDevTools
using Test

# write your own tests here
@time @testset "Benchmark Tests" begin
    include("benchmark_tests.jl")
end
@time @testset "ODE AppxTrue Tests" begin
    include("ode_appxtrue_tests.jl")
end
@time @testset "Analyticless Convergence Tests" begin
    include("analyticless_convergence_tests.jl")
end
@time @testset "ODE Tableau Convergence Tests" begin
    include("ode_tableau_convergence_tests.jl")
end ## Windows 32-bit fails on Butcher62 convergence test
@time @testset "Analyticless Stochastic WP" begin
    include("analyticless_stochastic_wp.jl")
end
@time @testset "Stability Region Tests" begin
    include("stability_region_test.jl")
end
@time @testset "Plot Recipes" begin
    include("plotrecipes_tests.jl")
end
@time @testset "Plot Recipes (Nonlinearsolve WP-diagrams)" begin
    include("nonlinearsolve_wpdiagram_tests.jl")
end
