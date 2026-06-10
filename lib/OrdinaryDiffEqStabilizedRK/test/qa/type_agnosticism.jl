using OrdinaryDiffEqStabilizedRK
using SciMLBase
using Test

# Checks that the stabilized RK methods do not promote Float32 problems to
# Float64 through hard-coded Float64 constants. Promotion of internal stage
# values is not observable from the solution (the typed integrator fields
# convert back), so instead `f` itself asserts the types it gets called with:
# any Float64 leaking into a stage value or stage time shows up here.
@testset "Float32 type agnosticism" begin
    all_solvers = [
        ROCK2(), ROCK4(), RKC(), RKMC2(), ESERK4(), ESERK5(), SERK2(),
        TSRKC2(), TSRKC3(), RKL1(), RKL2(), RKG1(), RKG2(),
    ]

    for solver in all_solvers
        @testset "$(nameof(typeof(solver))) out-of-place" begin
            u_ok = Ref(true)
            t_ok = Ref(true)
            f = (u, p, t) -> begin
                u isa Float32 || (u_ok[] = false)
                t isa Float32 || (t_ok[] = false)
                -20.0f0 * u
            end
            prob = ODEProblem(f, 1.0f0, (0.0f0, 1.0f0))
            sol = solve(prob, solver, abstol = 1.0f-4, reltol = 1.0f-4)
            @test SciMLBase.successful_retcode(sol)
            @test eltype(sol[end]) == Float32
            @test u_ok[]
            @test t_ok[]
        end

        @testset "$(nameof(typeof(solver))) in-place" begin
            u_ok = Ref(true)
            t_ok = Ref(true)
            f! = (du, u, p, t) -> begin
                eltype(u) == Float32 || (u_ok[] = false)
                t isa Float32 || (t_ok[] = false)
                du[1] = -20.0f0 * u[1]
                du[2] = -10.0f0 * u[2]
                nothing
            end
            prob = ODEProblem(f!, Float32[1.0, 1.0], (0.0f0, 1.0f0))
            sol = solve(prob, solver, abstol = 1.0f-4, reltol = 1.0f-4)
            @test SciMLBase.successful_retcode(sol)
            @test eltype(sol[end]) == Float32
            @test u_ok[]
            @test t_ok[]
        end
    end
end
