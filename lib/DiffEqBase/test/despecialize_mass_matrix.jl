using DiffEqBase, SciMLBase, LinearAlgebra, Test
using DiffEqBase: FunctionWrappersWrappers

# The resolved default algorithm uses finite differences, so concretization takes the
# non-ForwardDiff `promote_f` path. A mass matrix must not disable wrapping there,
# otherwise every mass-matrix model solved with the default algorithm keeps its own
# RHS type in the integrator and recompiles the whole solve.
struct FiniteDiffAlg <: SciMLBase.AbstractODEAlgorithm end
SciMLBase.forwarddiffs_model(::FiniteDiffAlg) = false
SciMLBase.isadaptive(::FiniteDiffAlg) = true

g1(du, u, p, t) = (du[1] = u[1]; du[2] = -u[2] + u[1]; nothing)
g2(du, u, p, t) = (du[1] = u[1]; du[2] = -2u[2] + u[1]; nothing)
M = Diagonal([0.0, 1.0])
u0 = [0.0, 1.0]

@testset "mass-matrix problems are wrapped on the finite-difference path" begin
    for spec in (SciMLBase.AutoDespecialize, SciMLBase.AutoSpecialize)
        probs = map((g1, g2)) do g
            ODEProblem{true, spec}(
                ODEFunction{true, spec}(g; mass_matrix = M), u0, (0.0, 1.0), (; k = 1.0)
            )
        end
        concrete = map(p -> DiffEqBase.get_concrete_problem(p, true; alg = FiniteDiffAlg()), probs)
        for c in concrete
            @test c.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper
            @test c.f.mass_matrix === M
            du = similar(u0)
            c.f(du, u0, c.p, 0.0)
            @test du[1] == u0[1]
        end
        @test typeof(concrete[1]) === typeof(concrete[2])
    end
end

@testset "a user Jacobian still opts out of the finite-difference wrapping path" begin
    jac!(J, u, p, t) = (J .= [1.0 0.0; 1.0 -1.0]; nothing)
    prob = ODEProblem{true, SciMLBase.AutoSpecialize}(
        ODEFunction{true, SciMLBase.AutoSpecialize}(g1; mass_matrix = M, jac = jac!),
        u0, (0.0, 1.0)
    )
    c = DiffEqBase.get_concrete_problem(prob, true; alg = FiniteDiffAlg())
    @test !(c.f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)
end
