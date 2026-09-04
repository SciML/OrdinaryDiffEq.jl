function _precompile_scalar_drift(u, p, t)
    return 1.01u
end

function _precompile_scalar_diffusion(u, p, t)
    return 0.87u
end

function _precompile_vector_drift!(du, u, p, t)
    du .= 1.01 .* u
    return nothing
end

function _precompile_vector_diffusion!(du, u, p, t)
    du .= 0.87 .* u
    return nothing
end

PrecompileTools.@compile_workload begin
    scalar_problem = SDEProblem(
        _precompile_scalar_drift, _precompile_scalar_diffusion, 1.0, (0.0, 0.1)
    )
    vector_problem = SDEProblem(
        _precompile_vector_drift!, _precompile_vector_diffusion!, ones(2), (0.0, 0.1)
    )
    solve(scalar_problem; dt = 0.1)
    solve(vector_problem; dt = 0.1)
end
