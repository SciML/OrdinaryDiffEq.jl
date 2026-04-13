using StochasticDiffEq, Test, Random
import SDEProblemLibrary: prob_sde_linear, prob_sde_2Dlinear
Random.seed!(100)

vecarrzero(m::Integer, n) = map(t -> zeros(n), 1:m)

tt = 0:0.05:1
ntt = length(tt)
out_VF = zeros(ntt)                                     # Vector{Float64}
out_VVF_1 = vecarrzero(ntt, 1)                           # Vector{Vector{Float64}}
out_VVF_2 = vecarrzero(ntt, 2)                           # Vector{Vector{Float64}}
out_VMF = vecarrzero(ntt, size(prob_sde_2Dlinear.u0))   # Vector{Matrix{Float64}}

@testset verbose = true "SDESolution interpolation" begin
    sol_SDE = solve(prob_sde_linear, SRIW1(); dt = 1 // 2^4)
    sol_SDE_2D = solve(prob_sde_2Dlinear, SRIW1(); dt = 1 // 2^4)

    sol_SDE_interp = sol_SDE(tt)
    sol_SDE_2D_interp = sol_SDE_2D(tt)

    @testset "1D" begin
        @test_throws MethodError sol_SDE(out_VF, tt; idxs = 1:1)
        @test sol_SDE(out_VF, tt) isa Vector{Float64}
        @test sol_SDE(out_VVF_1, tt) isa Vector{Vector{Float64}}
        @test sol_SDE_interp.u ≈ out_VF
    end

    @testset "2D" begin
        @test_throws MethodError sol_SDE_2D(out_VF, tt; idxs = 3:3)
        @test sol_SDE_2D(out_VF, tt; idxs = 3) isa Vector{Float64}
        @test sol_SDE_2D(out_VVF_1, tt; idxs = 3) isa Vector{Vector{Float64}}
        @test sol_SDE_2D(out_VVF_1, tt; idxs = 3:3) isa Vector{Vector{Float64}}
        @test sol_SDE_2D(out_VVF_2, tt; idxs = 2:3) isa Vector{Vector{Float64}}
        @test sol_SDE_2D(out_VMF, tt) isa Vector{Matrix{Float64}}
        @test sol_SDE_2D_interp.u ≈ out_VMF
    end
end
