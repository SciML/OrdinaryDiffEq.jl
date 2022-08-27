using OrdinaryDiffEq, CUDA, Test
CUDA.allowscalar(false)
N = 256
# Define the initial condition as normal arrays
u0 = zeros(N, N, 3)
u0 .= 1.0
gu0 = CuArray(Float32.(u0))

# Define the discretized PDE as an ODE function
function f(du, u, p, t)
    @. du = u
end
prob = ODEProblem(f, u0, (0.0f0, 10.0f0))
prob2 = ODEProblem(f, gu0, (0.0f0, 10.0f0))

algs = [ORK256(), CarpenterKennedy2N54(), SHLDDRK64(), DGLDDRK73_C()]

for alg in algs
    # CPU warmup
    solve(prob, alg, save_everystep = false, save_start = false, dt = 0.01)
    solve(prob, alg, save_everystep = false, save_start = false, dt = 0.01)
    solve(prob, alg, save_everystep = false, save_start = false, dt = 0.01)
    println("CPU Times for $alg")
    @time sol = solve(prob, alg, save_everystep = false, save_start = false, dt = 0.01)

    # GPU warmup
    solve(prob2, alg, save_everystep = false, save_start = false, dt = 0.01)
    solve(prob2, alg, save_everystep = false, save_start = false, dt = 0.01)
    solve(prob2, alg, save_everystep = false, save_start = false, dt = 0.01)
    println("GPU Times for $alg")
    @time sol2 = solve(prob2, alg, save_everystep = false, save_start = false, dt = 0.01)

    @test sol[end] ≈ Array(sol2[end])
end
