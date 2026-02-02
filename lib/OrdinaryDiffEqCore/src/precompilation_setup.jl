function lorenz(du, u, p, t)
    du[1] = 10.0(u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    return du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

function lorenz_oop(u, p, t)
    return [10.0(u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2], u[1] * u[2] - (8 / 3) * u[3]]
end

PrecompileTools.@compile_workload begin
    ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0))
    ODEProblem(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0), Float64[])
    ODEProblem{true, SciMLBase.AutoSpecialize}(
        lorenz, [1.0; 0.0; 0.0],
        (0.0, 1.0)
    )
    ODEProblem{true, SciMLBase.AutoSpecialize}(
        lorenz, [1.0; 0.0; 0.0],
        (0.0, 1.0), Float64[]
    )
    ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
        lorenz, [1.0; 0.0; 0.0],
        (0.0, 1.0)
    )
    ODEProblem{true, SciMLBase.FunctionWrapperSpecialize}(
        lorenz, [1.0; 0.0; 0.0],
        (0.0, 1.0), Float64[]
    )
    ODEProblem{true, SciMLBase.NoSpecialize}(lorenz, [1.0; 0.0; 0.0], (0.0, 1.0))
    ODEProblem{true, SciMLBase.NoSpecialize}(
        lorenz, [1.0; 0.0; 0.0], (0.0, 1.0),
        Float64[]
    )

    lorenz([1.0; 0.0; 0.0], [1.0; 0.0; 0.0], SciMLBase.NullParameters(), 0.0)
    lorenz([1.0; 0.0; 0.0], [1.0; 0.0; 0.0], Float64[], 0.0)
    lorenz_oop([1.0; 0.0; 0.0], SciMLBase.NullParameters(), 0.0)
    lorenz_oop([1.0; 0.0; 0.0], Float64[], 0.0)
end
