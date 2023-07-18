using OrdinaryDiffEq, Test
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear

vecarrzero(m::Integer, n) = map(t -> zeros(n), 1:m)

algs_ODE = zip([Tsit5(), Tsit5(), ABM54(), AutoTsit5(Rosenbrock23())],
    [Dict(), Dict(:dense => false), Dict(:dt => 0.1), Dict()],
    ["(normal)", "(not dense)", "(fixed-step)", "(CompositeAlgorithm)"])

tt = 0:0.05:1
ntt = length(tt)
out_VF = zeros(ntt)                                     # Vector{Float64}
out_VVF_1 = vecarrzero(ntt, 1)                           # Vector{Vector{Float64}}
out_VVF_2 = vecarrzero(ntt, 2)                           # Vector{Vector{Float64}}
out_VMF = vecarrzero(ntt, size(prob_ode_2Dlinear.u0))   # Vector{Matrix{Float64}}

@testset verbose=true "ODESolution interpolation $str" for (alg, kwargs, str) in algs_ODE
    sol_ODE = solve(prob_ode_linear, alg; kwargs...)
    sol_ODE_2D = solve(prob_ode_2Dlinear, alg; kwargs...)

    sol_ODE_interp = sol_ODE(tt)
    sol_ODE_2D_interp = sol_ODE_2D(tt)

    @testset "1D" begin
        @test_throws MethodError sol_ODE(out_VF, tt; idxs = 1:1)
        @test sol_ODE(out_VF, tt) isa Vector{Float64}
        @test sol_ODE(out_VVF_1, tt) isa Vector{Vector{Float64}}
        @test sol_ODE_interp.u == out_VF
    end

    @testset "2D" begin
        @test_throws MethodError sol_ODE_2D(out_VF, tt; idxs = 3:3)
        @test sol_ODE_2D(out_VF, tt; idxs = 3) isa Vector{Float64}
        @test sol_ODE_2D(out_VVF_1, tt; idxs = 3) isa Vector{Vector{Float64}}
        @test sol_ODE_2D(out_VVF_1, tt; idxs = 3:3) isa Vector{Vector{Float64}}
        @test sol_ODE_2D(out_VVF_2, tt; idxs = 2:3) isa Vector{Vector{Float64}}
        @test sol_ODE_2D(out_VMF, tt) isa Vector{Matrix{Float64}}
        @test sol_ODE_2D_interp.u == out_VMF
    end
end
