function _precompile_sde_drift!(du, u, p, t)
    du .= 1.01 .* u
    return nothing
end

function _precompile_sde_diffusion!(du, u, p, t)
    du .= 0.87 .* u
    return nothing
end

function _precompile_parameterized_sde_drift!(du, u, p, t)
    du .= p.drift .* u
    return nothing
end

function _precompile_parameterized_sde_diffusion!(du, u, p, t)
    du .= p.diffusion .* u
    return nothing
end

const _precompile_sde_parameters = (drift = 1.01, diffusion = 0.87)

PrecompileTools.@compile_workload begin
    problem_list = []
    if Preferences.@load_preference("PrecompileDefaultSpecialize", true)
        push!(
            problem_list,
            SDEProblem(
                _precompile_sde_drift!, _precompile_sde_diffusion!, ones(2), (0.0, 0.1)
            )
        )
        push!(
            problem_list,
            SDEProblem(
                _precompile_sde_drift!, _precompile_sde_diffusion!, ones(2),
                (0.0, 0.1), Float64[]
            )
        )
    end
    if Preferences.@load_preference("PrecompileAutoSpecialize", false)
        push!(
            problem_list,
            SDEProblem{true, SciMLBase.AutoSpecialize}(
                _precompile_sde_drift!, _precompile_sde_diffusion!, ones(2), (0.0, 0.1)
            )
        )
        push!(
            problem_list,
            SDEProblem{true, SciMLBase.AutoSpecialize}(
                _precompile_sde_drift!, _precompile_sde_diffusion!, ones(2),
                (0.0, 0.1), Float64[]
            )
        )
    end

    parameter_generic_problem_list = []
    if Preferences.@load_preference("PrecompileAutoDespecialize", true)
        push!(
            parameter_generic_problem_list,
            SDEProblem{true, SciMLBase.AutoDespecialize}(
                _precompile_parameterized_sde_drift!,
                _precompile_parameterized_sde_diffusion!, ones(2), (0.0, 0.1),
                _precompile_sde_parameters
            )
        )
    end
    if Preferences.@load_preference("PrecompileAutoDePSpecialize", true)
        push!(
            parameter_generic_problem_list,
            SDEProblem{true, SciMLBase.AutoDePSpecialize}(
                _precompile_parameterized_sde_drift!,
                _precompile_parameterized_sde_diffusion!, ones(2), (0.0, 0.1),
                _precompile_sde_parameters
            )
        )
    end
    append!(problem_list, parameter_generic_problem_list)

    if Preferences.@load_preference("PrecompileFunctionWrapperSpecialize", false)
        push!(
            problem_list,
            SDEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
                _precompile_sde_drift!, _precompile_sde_diffusion!, ones(2), (0.0, 0.1)
            )
        )
        push!(
            problem_list,
            SDEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
                _precompile_sde_drift!, _precompile_sde_diffusion!, ones(2),
                (0.0, 0.1), Float64[]
            )
        )
    end

    if Preferences.@load_preference("PrecompileNoSpecialize", false)
        push!(
            problem_list,
            SDEProblem{true, SciMLBase.NoSpecialize}(
                _precompile_sde_drift!, _precompile_sde_diffusion!, ones(2), (0.0, 0.1)
            )
        )
        push!(
            problem_list,
            SDEProblem{true, SciMLBase.NoSpecialize}(
                _precompile_sde_drift!, _precompile_sde_diffusion!, ones(2),
                (0.0, 0.1), Float64[]
            )
        )
    end

    for problem in problem_list
        solve(problem, SOSRI(); dt = 0.1)
    end
    for problem in parameter_generic_problem_list
        solve(problem; dt = 0.1)
    end
end
