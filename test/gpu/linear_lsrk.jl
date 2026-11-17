using OrdinaryDiffEq, CUDA, Test
using OrdinaryDiffEqLowStorageRK: ORK256, CarpenterKennedy2N54, SHLDDRK64, DGLDDRK73_C,
    CKLLSRK43_2, CKLLSRK54_3C, RDPK3SpFSAL35, RDPK3SpFSAL49
CUDA.allowscalar(false)
N = 256
# Define the initial condition as normal arrays
u0 = zeros(N, N, 3)
u0 .= 1.0
gu0 = CuArray(Float32.(u0))

# Define the discretized PDE as an ODE function
function f(du, u, p, t)
    return @. du = u
end
prob = ODEProblem{true, SciMLBase.FullSpecialize}(f, u0, (0.0f0, 10.0f0))
prob2 = ODEProblem{true, SciMLBase.FullSpecialize}(f, gu0, (0.0f0, 10.0f0))

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

    @test sol.u[end] ≈ Array(sol2.u[end])
end

@testset "Retried low-storage steps on the GPU" begin
    lin!(du, u, A, t) = (du .= A * u; nothing)
    prob = ODEProblem(lin!, CuArray([1.0, 0.0]), (0.0, 6.0), CuArray([0.0 1.0; -1.0 0.0]))
    algs = (CKLLSRK43_2(), CKLLSRK54_3C(), RDPK3SpFSAL35(), RDPK3SpFSAL49())
    for alg in algs, save_everystep in (true, false)
        integ = init(prob, alg; dt = 3.0, abstol = 1.0e-7, reltol = 1.0e-7, save_everystep)
        step!(integ)
        @test integ.stats.nreject > 0
        @test Array(integ.u) ≈ [cos(integ.t), -sin(integ.t)] atol = 2.0e-6
    end
end
