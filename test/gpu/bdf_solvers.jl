using OrdinaryDiffEqBDF, CUDA, Test
CUDA.allowscalar(false)

function f_gpu!(du, u, p, t)
    @. du = p * u
    return nothing
end

u0_cpu = [1.0, 2.0]
p_cpu = [-0.5, -1.5]
u0_gpu = cu(u0_cpu)
p_gpu = cu(p_cpu)
tspan = (0.0, 1.0)

prob_cpu = ODEProblem(f_gpu!, copy(u0_cpu), tspan, p_cpu)
prob_gpu = ODEProblem(f_gpu!, copy(u0_gpu), tspan, p_gpu)

@testset "QNDF GPU" begin
    sol_cpu = solve(prob_cpu, QNDF(), abstol = 1.0e-8, reltol = 1.0e-8)
    sol_gpu = solve(prob_gpu, QNDF(), abstol = 1.0e-8, reltol = 1.0e-8)
    @test sol_gpu.u[end] isa CuArray
    @test isapprox(Array(sol_gpu.u[end]), sol_cpu.u[end], rtol = 1.0e-5)
end

# FBDF GPU requires fixing scalar indexing in bdf_interpolants.jl (_get_theta)
@testset "FBDF GPU" begin
    @test_broken begin
        sol_gpu = solve(prob_gpu, FBDF(), abstol = 1.0e-8, reltol = 1.0e-8)
        sol_cpu = solve(prob_cpu, FBDF(), abstol = 1.0e-8, reltol = 1.0e-8)
        isapprox(Array(sol_gpu.u[end]), sol_cpu.u[end], rtol = 1.0e-5)
    end
end
