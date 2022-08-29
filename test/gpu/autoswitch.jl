using OrdinaryDiffEq, CUDA, Test
CUDA.allowscalar(false)

# https://github.com/SciML/OrdinaryDiffEq.jl/issues/1614
function f(du, u, p, t)
    @. du = u
end

problem = ODEProblem(f, CUDA.ones(1), (0.0, 1.0))
stiffalg = Rosenbrock23()
for alg in (AutoDP5, AutoTsit5, AutoVern6, AutoVern7, AutoVern8, AutoVern9)
    @test_nowarn solve(problem, alg(stiffalg), save_everystep=false, save_start=false)
end

