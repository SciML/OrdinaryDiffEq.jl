using OrdinaryDiffEqLowOrderRK: Euler
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqVerner: Vern9
using DiffEqBase, DiffEqDevTools, Test
using ODEProblemLibrary: prob_ode_linear

f = (u, p, t) -> u
fun = ODEFunction(f; analytic = (u0, p, t) -> u0 * exp(t))
prob = ODEProblem(fun, 1 / 2, (0.0, 1.0))

sol = solve(prob, Euler(); dt = 1 // 2^(4), dense_errors = true)
sol2 = solve(prob, Vern9(); dt = 1 // 2^(10), abstol = 1.0e-14, reltol = 1.0e-14)

prob2 = ODEProblem(f, 1 / 2, (0.0, 1.0))
sol3 = solve(prob_ode_linear, Euler(); dt = 1 // 2^(4))

errsol1 = appxtrue(sol, sol2)

sol4 = solve(prob, Euler(); dt = 1 // 2^(4))
test_sol = TestSolution(sol2)
errsol2 = appxtrue(sol4, test_sol)

sol5 = solve(prob, Euler(); dt = 1 // 2^(4))
test_sol = TestSolution(sol2.t, sol2.u[end])
errsol3 = appxtrue(sol5, test_sol)

@test errsol1.errors[:L2] ≈ 0.018865798306718855
@test errsol1.errors[:L2] ≈ errsol2.errors[:L2]

@test errsol1.errors[:final] ≈ sol.errors[:final]
@test errsol1.errors[:l2] ≈ sol.errors[:l2]
@test errsol1.errors[:l∞] ≈ sol.errors[:l∞]
@test errsol1.errors[:L∞] ≈ sol.errors[:L∞]
@test errsol1.errors[:L2] ≈ sol.errors[:L2]

sol_nodense = solve(
    prob, Vern9(); dt = 1 // 2^(10), abstol = 1.0e-14, reltol = 1.0e-14,
    dense = false, saveat = 1 // 4
)

@test TestSolution(sol2).dense
@test !TestSolution(sol_nodense).dense
@test haskey(appxtrue(sol4, TestSolution(sol_nodense)).errors, :l2)

@testset "Nondense TestSolution sampled errors" begin
    for u0 in (1.0, [1.0, 2.0], [1.0 2.0; 3.0 4.0])
        @testset "$(typeof(u0))" begin
            problem = ODEProblem((u, p, t) -> -u, u0, (0.0, 1.0))
            reference = solve(problem, Tsit5(); abstol = 1.0e-12, reltol = 1.0e-12)
            approximate = solve(problem, Tsit5(); dense = false, adaptive = false, dt = 0.25)
            wrapper = TestSolution(reference)
            wrapper.dense = false
            result = appxtrue(approximate, wrapper)
            sampled = reference(approximate.t).u
            differences = [a - b for (a, b) in zip(approximate.u, sampled)]
            flat = u0 isa Number ? differences : reduce(vcat, vec.(differences))
            @test Set(keys(result.errors)) == Set((:final, :l2, :l∞))
            @test result.errors[:l∞] ≈ maximum(abs, flat)
            @test result.errors[:l2] ≈ sqrt(sum(abs2, flat) / length(flat))
            @test result.errors[:final] ≈ sum(abs, differences[end]) / length(differences[end])

            # A reference constructed without dense output is still sampled on the query grid.
            sparse_reference = solve(problem, Tsit5(); dense = false, saveat = 0.5)
            sparse_wrapper = TestSolution(sparse_reference)
            @test !sparse_wrapper.dense
            sparse_result = appxtrue(approximate, sparse_wrapper)
            sparse_samples = sparse_reference(approximate.t).u
            sparse_differences = [a - b for (a, b) in zip(approximate.u, sparse_samples)]
            sparse_flat = u0 isa Number ? sparse_differences : reduce(vcat, vec.(sparse_differences))
            @test Set(keys(sparse_result.errors)) == Set((:final, :l2, :l∞))
            @test sparse_result.errors[:l∞] ≈ maximum(abs, sparse_flat)
            @test sparse_result.errors[:l2] ≈ sqrt(sum(abs2, sparse_flat) / length(sparse_flat))

            stored = TestSolution(copy(approximate.t), deepcopy(approximate.u))
            @test all(iszero, values(appxtrue(approximate, stored).errors))
            mismatched = TestSolution(approximate.t[1:(end - 1)], approximate.u[1:(end - 1)])
            @test Set(keys(appxtrue(approximate, mismatched).errors)) == Set((:final,))
            dense_result = appxtrue(approximate, TestSolution(reference))
            @test Set(keys(dense_result.errors)) == Set((:final, :l2, :l∞, :L2, :L∞))
            @test dense_result.errors[:l2] ≈ result.errors[:l2]
            @test dense_result.errors[:l∞] ≈ result.errors[:l∞]
        end
    end
end
