using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using OrdinaryDiffEqNonlinearSolve: BrownFullBasicInit
using CUDA
using LinearAlgebra
using Adapt
using SparseArrays
using Test

#=
du[1] = -u[1]
du[2] = -0.5*u[2]
    0 =  u[1] + u[2] - u[3]
    0 = -u[1] + u[2] - u[4]
=#

function dae!(du, u, p, t)
    return mul!(du, p, u)
end

p = [
    -1 0 0 0
    1 -0.5 0 0
    1 1 -1 0
    -1 1 0 -1
]

# mass_matrix = [1 0 0 0
#                0 1 0 0
#                0 0 0 0
#                0 0 0 0]
mass_matrix = Diagonal([1, 1, 0, 0])
jac_prototype = sparse(map(x -> iszero(x) ? 0.0 : 1.0, p))

u0 = [1.0, 1.0, 0.5, 0.5] # force init
odef = ODEFunction(dae!, mass_matrix = mass_matrix, jac_prototype = jac_prototype)

tspan = (0.0, 5.0)
prob = ODEProblem(odef, u0, tspan, p)
# v7: DefaultInit is CheckInit, so pass BrownFullBasicInit explicitly to let
# the solver fix the inconsistent initial u0 before integrating.
sol = solve(prob, Rodas5P(); initializealg = BrownFullBasicInit())

# gpu version
mass_matrix_d = adapt(CuArray, mass_matrix)

# TODO: jac_prototype fails
# jac_prototype_d = adapt(CuArray, jac_prototype)
# jac_prototype_d = CUDA.CUSPARSE.CuSparseMatrixCSR(jac_prototype)
jac_prototype_d = nothing

u0_d = adapt(CuArray, u0)
p_d = adapt(CuArray, p)
odef_d = ODEFunction(dae!, mass_matrix = mass_matrix_d, jac_prototype = jac_prototype_d)
prob_d = ODEProblem(odef_d, u0_d, tspan, p_d)
sol_d = solve(prob_d, Rodas5P(); initializealg = BrownFullBasicInit())

@testset "Test constraints in GPU sol" begin
    for t in sol_d.t
        u = Vector(sol_d(t))
        @test isapprox(u[1] + u[2], u[3]; atol = 1.0e-6)
        @test isapprox(-u[1] + u[2], u[4]; atol = 1.0e-6)
    end
end

@testset "Compare GPU to CPU solution" begin
    for t in tspan[begin]:0.1:tspan[end]
        @test isapprox(Vector(sol_d(t)), sol(t); rtol = 1.0e-4)
    end
end
