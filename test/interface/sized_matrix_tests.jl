using OrdinaryDiffEq, StaticArrays, Test

function dudt!(du, σ, p, t)
    return du .= -σ
end

η0 = 1
τ = 1
α = 0.8
p_giesekus = [η0, τ, α]

## Exercise code paths that are in-place but not `<:Array` but not AbstractVector
σ0 = SizedMatrix{3, 3}([1.0 0.0 0.0; 0.0 2.0 0.0; 3.0 0.0 0.0])

prob_giesekus = ODEProblem(dudt!, σ0, (0.0, 2.0), p_giesekus)

solve_giesekus = solve(
    prob_giesekus, Rodas4(), saveat = 0.2, abstol = 1.0e-14,
    reltol = 1.0e-14
)
for alg in [
        Rosenbrock23(),
        Rodas4(),
        Rodas4P(),
        Rodas5(),
        Rodas5P(),
        Tsit5(),
        Vern6(),
        Vern7(),
        Vern8(),
        Vern9(),
        DP5(),
    ]
    sol = solve(prob_giesekus, alg, saveat = 0.2, abstol = 1.0e-14, reltol = 1.0e-14)
    @test Array(sol) ≈ Array(solve_giesekus)
end
