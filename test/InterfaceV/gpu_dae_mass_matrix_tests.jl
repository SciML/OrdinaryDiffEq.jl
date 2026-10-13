using OrdinaryDiffEq, OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK,
    OrdinaryDiffEqBDF, OrdinaryDiffEqNonlinearSolve,
    JLArrays, LinearAlgebra, SparseArrays, Test
import DiffEqBase

#=
Regression test for the Julia-1.13 GPU DAE regression: `==(A::AbstractMatrix,
::UniformScaling)` reads `first(A)` on 1.13+, so `f.mass_matrix != I` in
`_ode_init` (and `all(isequal(0), f.mass_matrix)` for non-diagonal-only
Rosenbrock methods) scalar-indexed a GPU-backed mass matrix. `JLArray` emulates
GPU scalar-indexing semantics on CPU, so a `Diagonal{T, <:JLArray}` mass matrix
exercises the same code path as `cu(Diagonal(...))` on CUDA.

Unlike CUDA, JLArrays does not convert `AbstractArray{Bool}` indices to integer
indices in `to_index`/`to_indices` (CUDA.jl does `to_index(::CuArray, I) =
findall(I)`). Emulate that rule here so `view(u, bool_mask)` in
`BrownFullBasicInit` works the same way it does on CUDA.
=#
Base.to_index(::JLArray, I::AbstractArray{Bool}) = JLArray(findall(Array(I)))
if VERSION >= v"1.11.0-DEV.1157"
    Base.to_indices(A::JLArray, I::Tuple{AbstractArray{Bool}}) =
        (Base.to_index(A, I[1]),)
end

@testset "GPU-emulated (JLArray) diagonal-mass-matrix DAE" begin
    function dae!(du, u, p, t)
        return mul!(du, p, u)
    end

    P = [
        -1 0 0 0
        1 -0.5 0 0
        1 1 -1 0
        -1 1 0 -1
    ]
    MASS_MATRIX = Diagonal([1, 1, 0, 0])
    JAC_PROTOTYPE = sparse(map(x -> iszero(x) ? 0.0 : 1.0, P))
    U0 = [1.0, 1.0, 0.5, 0.5]
    TSPAN = (0.0, 5.0)

    odef_cpu = ODEFunction{true, SciMLBase.FullSpecialize}(
        dae!; mass_matrix = MASS_MATRIX,
        jac_prototype = JAC_PROTOTYPE
    )
    prob_cpu = ODEProblem(
        odef_cpu, U0, TSPAN, P;
        initializealg = DiffEqBase.BrownFullBasicInit()
    )

    M_D = jl(MASS_MATRIX)
    odef_d = ODEFunction{true, SciMLBase.FullSpecialize}(
        dae!; mass_matrix = M_D,
        jac_prototype = nothing
    )
    prob_d = ODEProblem(
        odef_d, jl(U0), TSPAN, jl(P);
        initializealg = DiffEqBase.BrownFullBasicInit()
    )

    for ALG in (
            Rosenbrock23, ROS2, Rodas5P, SDIRK2, Cash4, Hairer4,
            ABDF2, QNDF2, QBDF2, FBDF, ImplicitEuler,
        )
        sol_cpu = solve(prob_cpu, ALG(); maxiters = 10_000)
        sol_d = solve(prob_d, ALG(); maxiters = 10_000)
        @test sol_d.retcode == SciMLBase.ReturnCode.Success
        @test Array(sol_d.u[end]) ≈ sol_cpu.u[end] atol = 1.0e-3
    end
end
