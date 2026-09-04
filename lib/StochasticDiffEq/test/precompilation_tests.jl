using StochasticDiffEq, Test

scalar_drift(u, p, t) = 1.01u
scalar_diffusion(u, p, t) = 0.87u

function vector_drift!(du, u, p, t)
    du .= 1.01 .* u
    return nothing
end

function vector_diffusion!(du, u, p, t)
    du .= 0.87 .* u
    return nothing
end

scalar_problem = SDEProblem(
    scalar_drift, scalar_diffusion, 1.0, (0.0, 0.1)
)
vector_problem = SDEProblem(
    vector_drift!, vector_diffusion!, ones(2), (0.0, 0.1)
)

@test solve(scalar_problem; dt = 0.1).retcode == ReturnCode.Success
@test solve(vector_problem; dt = 0.1).retcode == ReturnCode.Success
@test init(scalar_problem; dt = 0.1).alg isa SOSRI
@test init(vector_problem; dt = 0.1).alg isa SOSRI
