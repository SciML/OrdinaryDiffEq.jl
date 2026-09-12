using OrdinaryDiffEqCore, OrdinaryDiffEqTsit5, SciMLBase, LinearAlgebra, Enzyme, Test

const TEST_A = [-0.3 1.0; -1.0 -0.3]
function test_base!(du, u, p, t)
    mul!(du, TEST_A, u)
    return nothing
end
const TEST_REF = solve(
    ODEProblem(test_base!, [1.0, 0.5], (0.0, 3.0)), Tsit5();
    abstol = 1.0e-10, reltol = 1.0e-10, dense = true, save_everystep = true
)

function test_forced!(du, u, p, t)
    mul!(du, TEST_A, u)
    r = TEST_REF(t)
    du[1] += p[1] * r[1]
    du[2] += p[1] * r[2]
    return nothing
end

@testset "Enzyme reverse through interpolant" begin
    u = [1.0, 0.5]
    p = [0.0]
    t = 1.5
    du = zeros(2)
    bdu = [1.0, 0.0]
    bu = zeros(2)
    bp = zeros(1)
    Enzyme.autodiff(
        Enzyme.Reverse, test_forced!, Enzyme.Const,
        Enzyme.Duplicated(du, copy(bdu)),
        Enzyme.Duplicated(u, bu),
        Enzyme.Duplicated(p, bp),
        Enzyme.Const(t)
    )
    @test bu ≈ [-0.3, 1.0]
    @test bp[1] ≈ TEST_REF(t)[1]
end
