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

function other_vector_drift!(du, u, p, t)
    du .= 0.99 .* u
    return nothing
end

function other_vector_diffusion!(du, u, p, t)
    du .= 0.76 .* u
    return nothing
end

function parameterized_vector_drift!(du, u, p, t)
    du .= p.drift .* u
    return nothing
end

function parameterized_vector_diffusion!(du, u, p, t)
    du .= p.diffusion .* u
    return nothing
end

function other_parameterized_vector_drift!(du, u, p, t)
    du .= (p.drift - 0.01) .* u
    return nothing
end

function other_parameterized_vector_diffusion!(du, u, p, t)
    du .= (p.diffusion - 0.01) .* u
    return nothing
end

function matrix_diffusion!(du, u, p, t)
    du[:, 1] .= 0.87 .* u
    return nothing
end

scalar_problem = SDEProblem(
    scalar_drift, scalar_diffusion, 1.0, (0.0, 0.1)
)
vector_problem = SDEProblem(
    vector_drift!, vector_diffusion!, ones(2), (0.0, 0.1)
)
other_vector_problem = SDEProblem(
    other_vector_drift!, other_vector_diffusion!, ones(2), (0.0, 0.1)
)

@test solve(scalar_problem; dt = 0.1).retcode == ReturnCode.Success
@test solve(vector_problem; dt = 0.1).retcode == ReturnCode.Success
@test init(scalar_problem; dt = 0.1).alg isa SOSRI
@test init(vector_problem; dt = 0.1).alg isa SOSRI

vector_f = init(vector_problem, SOSRI(); dt = 0.1).sol.prob.f
other_vector_f = init(other_vector_problem, SOSRI(); dt = 0.1).sol.prob.f
@test SciMLBase.specialization(vector_f) === SciMLBase.AutoSpecialize
@test typeof(vector_f) === typeof(other_vector_f)

parameters = (drift = 1.01, diffusion = 0.87)
for specialization in (
        SciMLBase.AutoDespecialize,
        SciMLBase.AutoRespecialize,
        SciMLBase.FunctionWrapperSpecialize,
        SciMLBase.NoSpecialize,
    )
    parameterized_problem = SDEProblem{true, specialization}(
        parameterized_vector_drift!, parameterized_vector_diffusion!, ones(2),
        (0.0, 0.1), parameters
    )
    other_parameterized_problem = SDEProblem{true, specialization}(
        other_parameterized_vector_drift!, other_parameterized_vector_diffusion!, ones(2),
        (0.0, 0.1), parameters
    )
    parameterized_f = init(parameterized_problem, SOSRI(); dt = 0.1).sol.prob.f
    other_parameterized_f = init(other_parameterized_problem, SOSRI(); dt = 0.1).sol.prob.f
    @test SciMLBase.specialization(parameterized_f) === specialization
    @test typeof(parameterized_f) === typeof(other_parameterized_f)
    @test solve(parameterized_problem; dt = 0.1).retcode == ReturnCode.Success
end

matrix_problem = SDEProblem{true, SciMLBase.AutoSpecialize}(
    vector_drift!, matrix_diffusion!, ones(2), (0.0, 0.1);
    noise_rate_prototype = zeros(2, 1)
)
@test solve(matrix_problem, EM(); dt = 0.01).retcode == ReturnCode.Success
@test solve(vector_problem, ImplicitRKMil(); dt = 0.01, adaptive = false).retcode ==
    ReturnCode.Success
