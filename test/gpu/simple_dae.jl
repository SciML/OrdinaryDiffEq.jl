using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using CUDA
using LinearAlgebra
using Adapt
using SparseArrays
using Test
using CUDSS

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
sol = solve(prob, Rodas5P())

mass_matrix_d = adapt(CuArray, mass_matrix)
jac_prototype_d_csc = CUDA.CUSPARSE.CuSparseMatrixCSC(jac_prototype)
jac_prototype_d_csr = CUDA.CUSPARSE.CuSparseMatrixCSR(jac_prototype)

u0_d = adapt(CuArray, u0)
p_d = adapt(CuArray, p)

odef_d_dns = ODEFunction(dae!, mass_matrix = mass_matrix_d, jac_prototype = nothing)
odef_d_csc = ODEFunction(dae!, mass_matrix = mass_matrix_d, jac_prototype = jac_prototype_d_csc)
odef_d_csr = ODEFunction(dae!, mass_matrix = mass_matrix_d, jac_prototype = jac_prototype_d_csr)

prob_d_dns = ODEProblem(odef_d_dns, u0_d, tspan, p_d)
prob_d_csc = ODEProblem(odef_d_csc, u0_d, tspan, p_d)
prob_d_csr = ODEProblem(odef_d_csr, u0_d, tspan, p_d)

sol_d_dns = solve(prob_d_dns, Rodas5P()) # works
sol_d_csc = solve(prob_d_csc, Rodas5P()) # uses Krylov
sol_d_csr = solve(prob_d_csr, Rodas5P()) # uses CUDSS

@testset "Test constraints in GPU sol" begin
    for sol_d in [sol_d_dns, sol_d_csc, sol_d_csr]
        for t in sol_d.t
            u = Vector(sol_d(t))
            @test isapprox(u[1] + u[2], u[3]; atol = 1.0e-6)
            @test isapprox(-u[1] + u[2], u[4]; atol = 1.0e-6)
        end
    end
end

@testset "Compare GPU to CPU solution" begin
    for sol_d in [sol_d_dns, sol_d_csc, sol_d_csr]
        for t in tspan[begin]:0.1:tspan[end]
            @test Vector(sol_d(t)) ≈ sol(t) atol=1e-5
        end
    end
end
