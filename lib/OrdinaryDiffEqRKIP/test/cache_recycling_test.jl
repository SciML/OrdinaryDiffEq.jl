using SciMLOperators: MatrixOperator
using OrdinaryDiffEqRKIP: RKIP
using SciMLBase: SplitODEProblem, SplitFunction, solve

import LinearAlgebra

@testset "Cache Recycling Test" begin
    matrix = [1 0.2 -0.2; 0.4 1 -1; 0.2 0.0 1.2]
    Â = MatrixOperator(matrix)
    u0 = [0.0, 0.0, 0.0]
    f! = (dx, x, p, t) -> dx .= cos.(x)

    alg = RKIP(1e-2, 1e1)

    splfc = SplitFunction{true}(Â, f!)
    spltode = SplitODEProblem(splfc, u0, (0, 1.0))

    sol_1 = solve(spltode, alg).u[end]
    sol_2 = solve(spltode, alg).u[end]
    @test all(sol_1 .≈ sol_2)
end
