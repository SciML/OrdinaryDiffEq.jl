using SafeTestsets

@time @safetestset "Linear-Nonlinear Krylov Methods Tests" include("linear_nonlinear_krylov_tests.jl")
@time @safetestset "Linear-Nonlinear Convergence Tests" include("linear_nonlinear_convergence_tests.jl")
