using OrdinaryDiffEqBDF, CUDA, Test
CUDA.allowscalar(false)

function f_gpu!(du, u, p, t)
    @. du = p * u
    return nothing
end

u0_cpu = [1.0, 2.0]
p_cpu = [-0.5, -1.5]
u0_gpu = CuArray(u0_cpu)
p_gpu = CuArray(p_cpu)
tspan = (0.0, 1.0)

prob_cpu = ODEProblem(f_gpu!, copy(u0_cpu), tspan, p_cpu)
prob_gpu = ODEProblem(f_gpu!, copy(u0_gpu), tspan, p_gpu)

@testset "QNDF GPU" begin
    sol_cpu = solve(prob_cpu, QNDF(), abstol = 1.0e-6, reltol = 1.0e-6)
    sol_gpu = solve(prob_gpu, QNDF(), abstol = 1.0e-6, reltol = 1.0e-6)
    @test sol_gpu.u[end] isa CuArray
    @test isapprox(Array(sol_gpu.u[end]), sol_cpu.u[end], rtol = 1.0e-3)
end

# FBDF GPU - scalar indexing in _bdf_active_order fixed via _is_all_zero helper
@testset "FBDF GPU" begin
    sol_cpu = solve(prob_cpu, FBDF(), abstol = 1.0e-6, reltol = 1.0e-6)
    sol_gpu = solve(prob_gpu, FBDF(), abstol = 1.0e-6, reltol = 1.0e-6)
    @test sol_gpu.u[end] isa CuArray
    @test isapprox(Array(sol_gpu.u[end]), sol_cpu.u[end], rtol = 1.0e-3)
end
