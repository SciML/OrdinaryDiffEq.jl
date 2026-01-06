using SafeTestsets

@time @safetestset "Low Order ERK Convergence Tests" include("low_order_erk_convergence_tests.jl")
@time @safetestset "OwrenZen Tests" include("owrenzen_tests.jl")
@time @safetestset "Euler SSP Tests" include("euler_ssp.jl")

# Only run QA and allocation tests on stable Julia versions
if isempty(VERSION.prerelease)
    @time @safetestset "JET Tests" include("jet.jl")
    @time @safetestset "Aqua" include("qa.jl")
    @time @safetestset "Allocation Tests" include("allocation_tests.jl")
end
