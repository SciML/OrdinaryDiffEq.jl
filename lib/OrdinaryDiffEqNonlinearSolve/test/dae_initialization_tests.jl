using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, OrdinaryDiffEqNonlinearSolve,
    StaticArrays, LinearAlgebra, Test, ADTypes, ForwardDiff, SciMLBase
using OrdinaryDiffEqNonlinearSolve: default_nlsolve

## Mass Matrix

function rober_oop(u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du1 = -k₁ * y₁ + k₃ * y₂ * y₃
    du2 = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du3 = y₁ + y₂ + y₃ - 1
    return [du1, du2, du3]
end

M = Diagonal([1.0, 1.0, 0.0])
f_oop = ODEFunction(rober_oop, mass_matrix = M)
prob_mm = ODEProblem(f_oop, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
sol = @inferred solve(
    prob_mm, Rosenbrock23(autodiff = AutoFiniteDiff()), reltol = 1.0e-8, abstol = 1.0e-8
)
@test sol.u[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!
sol = @inferred solve(
    prob_mm, Rosenbrock23(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = ShampineCollocationInit()
)
@test sol.u[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!

integrator = @inferred init(prob_mm, Rodas5(autodiff = AutoForwardDiff(chunksize = 3)))
# It would be nice if this could test that the initialization is fully inferred,
# but since the return is just the integrator, this doesn't meet that goal.
@inferred SciMLBase.initialize_dae!(integrator)

prob_mm = ODEProblem(f_oop, [1.0, 0.0, 0.2], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
sol = solve(
    prob_mm, Rosenbrock23(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)
@test sum(sol.u[1]) ≈ 1
@test sol.u[1] ≈ [1.0, 0.0, 0.0]
for alg in [Rosenbrock23(autodiff = AutoFiniteDiff()), Trapezoid()]
    local sol
    sol = solve(
        prob_mm, alg, reltol = 1.0e-8, abstol = 1.0e-8,
        initializealg = ShampineCollocationInit()
    )
    @test sum(sol.u[1]) ≈ 1
end

function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end
M = Diagonal([1.0, 1.0, 0.0])
f = ODEFunction(rober, mass_matrix = M)
prob_mm = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
sol = @inferred solve(prob_mm, Rodas5(autodiff = AutoFiniteDiff()), reltol = 1.0e-8, abstol = 1.0e-8)
@test sol.u[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!
sol = solve(
    prob_mm, Rodas5(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = ShampineCollocationInit()
)
@test sol.u[1] == [1.0, 0.0, 0.0] # Ensure initialization is unchanged if it works at the start!

integrator = @inferred init(prob_mm, Rodas5(autodiff = AutoForwardDiff(chunksize = 3)))
# It would be nice if this could test that the initialization is fully inferred,
# but since the return is just the integrator, this doesn't meet that goal.
@inferred SciMLBase.initialize_dae!(integrator)

prob_mm = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1.0e5), (0.04, 3.0e7, 1.0e4))
sol = solve(
    prob_mm, Rodas5(), reltol = 1.0e-8, abstol = 1.0e-8,
    initializealg = BrownFullBasicInit()
)
@test sum(sol.u[1]) ≈ 1
@test sol.u[1] ≈ [1.0, 0.0, 0.0]

for alg in [Rodas5(autodiff = AutoFiniteDiff()), Trapezoid()]
    local sol
    sol = solve(
        prob_mm, alg, reltol = 1.0e-8, abstol = 1.0e-8,
        initializealg = ShampineCollocationInit()
    )
    @test sum(sol.u[1]) ≈ 1
end

function rober_no_p(du, u, p, t)
    y₁, y₂, y₃ = u
    (k₁, k₂, k₃) = (0.04, 3.0e7, 1.0e4)
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    return nothing
end

function rober_oop_no_p(du, u, p, t)
    y₁, y₂, y₃ = u
    (k₁, k₂, k₃) = (0.04, 3.0e7, 1.0e4)
    du1 = -k₁ * y₁ + k₃ * y₂ * y₃
    du2 = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du3 = y₁ + y₂ + y₃ - 1
    return [du1, du2, du3]
end

# test oop and iip ODE initialization with parameters without eltype/length
struct UnusedParam
end
for f in (
        ODEFunction(rober_no_p, mass_matrix = M), ODEFunction(rober_oop_no_p, mass_matrix = M),
    )
    local prob, probp
    prob = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1.0e5))
    probp = ODEProblem(f, [1.0, 0.0, 1.0], (0.0, 1.0e5), UnusedParam)
    for initializealg in (ShampineCollocationInit(), BrownFullBasicInit())
        isapprox(
            init(prob, Rodas5(), abstol = 1.0e-10; initializealg).u,
            init(probp, Rodas5(), abstol = 1.0e-10; initializealg).u
        )
    end
end

# to test that we get the right NL solve we need a broken solver.
struct BrokenNLSolve <: SciMLBase.AbstractNonlinearAlgorithm
    BrokenNLSolve(; kwargs...) = new()
end
function SciMLBase.__solve(
        prob::NonlinearProblem,
        alg::BrokenNLSolve, args...;
        kwargs...
    )
    u = fill(reinterpret(Float64, 0xDEADBEEFDEADBEEF), 3)
    return SciMLBase.build_solution(
        prob, alg, u, copy(u);
        retcode = ReturnCode.Success
    )
end
function f2(u, p, t)
    return u
end
f = ODEFunction(f2, mass_matrix = Diagonal([1.0, 1.0, 0.0]))
prob = ODEProblem(f, ones(3), (0.0, 1.0))
integrator = init(
    prob, Rodas5P(),
    initializealg = ShampineCollocationInit(1.0, BrokenNLSolve())
)
@test all(isequal(reinterpret(Float64, 0xDEADBEEFDEADBEEF)), integrator.u)

@testset "`reinit!` reruns initialization" begin
    initializeprob = NonlinearProblem(1.0, [0.0]) do u, p
        return u^2 - p[1]^2
    end
    initializeprobmap = function (nlsol)
        return [nlsol.prob.p[1], nlsol.u]
    end
    update_initializeprob! = function (iprob, integ)
        iprob.p[1] = integ.u[1]
    end
    initialization_data = SciMLBase.OverrideInitData(
        initializeprob, update_initializeprob!, initializeprobmap, nothing
    )
    fn = ODEFunction(; mass_matrix = [1 0; 0 0], initialization_data) do du, u, p, t
        du[1] = u[1]
        du[2] = u[1]^2 - u[2]^2
    end
    prob = ODEProblem(fn, [2.0, 0.0], (0.0, 1.0))
    integ = init(prob, Rodas5P(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test integ.u ≈ [2.0, 2.0] atol = 1.0e-8
    reinit!(integ)
    @test integ.u ≈ [2.0, 2.0] atol = 1.0e-8
    step!(integ, 0.01, true)
    @test SciMLBase.successful_retcode(integ.sol.retcode)
    reinit!(integ, reinit_dae = false)
    @test integ.u ≈ [2.0, 0.0]
    # With reinit_dae=false the algebraic constraint u[1]^2 - u[2]^2 = 0 is violated.
    # Rosenbrock methods (Rodas5P) linearize and don't iterate, so the step succeeds
    # but the constraint remains violated — u[2] stays near 0 instead of tracking u[1].
    step!(integ, 0.01, true)
    @test abs(integ.u[2]) < 1.0e-10  # u[2] stuck near 0, not reinitialized
    @test abs(integ.u[1]) > 1.5    # u[1] still evolving
end

function override_dual_init_residual!(residual, z, parameters)
    c, x = parameters
    residual[1] = z[1]^3 + z[1] - c * x
    return nothing
end

function override_dual_init_residual(z, parameters)
    c, x = parameters
    return [z[1]^3 + z[1] - c * x]
end

function override_dual_ode_rhs!(du, u, p, t)
    c, k = p
    x, z = u
    du[1] = -k * x + z
    du[2] = z^3 + z - c * x
    return nothing
end

function override_dual_ode_rhs(u, p, t)
    c, k = p
    x, z = u
    return [-k * x + z, z^3 + z - c * x]
end

function make_override_dual_problem(p, inplace, vector_init_params)
    c, _ = p
    x0 = 0.5 * c
    init_params = vector_init_params ? [c, x0] : (c, x0)
    initprob = if inplace
        NonlinearProblem(override_dual_init_residual!, [one(c) * 0.5], init_params)
    else
        NonlinearProblem(override_dual_init_residual, [one(c) * 0.5], init_params)
    end
    update_initializeprob! = function (ip, valp)
        c2, _ = SciMLBase.parameter_values(valp)
        x2 = 0.5 * c2
        new_params = vector_init_params ? [c2, x2] : (c2, x2)
        return SciMLBase.remake(ip; u0 = [one(c2) * 0.5], p = new_params)
    end
    initializeprobmap = nlsol -> begin
        _, x2 = nlsol.prob.p
        [x2, nlsol.u[1]]
    end
    initialization_data = SciMLBase.OverrideInitData(
        initprob, update_initializeprob!, initializeprobmap, nothing;
        is_update_oop = Val(true)
    )
    mass_matrix = Diagonal([1.0, 0.0])
    f = if inplace
        ODEFunction{true}(
            override_dual_ode_rhs!; mass_matrix, initialization_data
        )
    else
        ODEFunction{false}(
            override_dual_ode_rhs; mass_matrix, initialization_data
        )
    end
    T = promote_type(eltype(p), Float64)
    return ODEProblem(f, T[ForwardDiff.value(x0), 0.5], (0.0, 1.0), p)
end

@testset "OverrideInit preserves algebraic ForwardDiff partials" begin
    tolerances = (; abstol = 1.0e-12, reltol = 1.0e-12)
    expected_state = [1.0, 1.0]
    expected_jacobian = [0.5 0.0; 0.5 0.0]
    for inplace in (true, false), vector_init_params in (false, true)
        label = "$(inplace ? "in-place" : "out-of-place"), " *
            "$(vector_init_params ? "vector" : "tuple") initialization parameters"
        @testset "$label" begin
            p = [2.0, 0.5]
            prob = make_override_dual_problem(p, inplace, vector_init_params)
            integrator = init(prob, Rodas5P(); tolerances...)
            @test integrator.u ≈ expected_state atol = 1.0e-8
            init_jacobian = ForwardDiff.jacobian(p) do dual_p
                dual_prob = make_override_dual_problem(
                    dual_p, inplace, vector_init_params
                )
                dual_integrator = init(dual_prob, Rodas5P(); tolerances...)
                collect(dual_integrator.u)
            end
            @test init_jacobian ≈ expected_jacobian atol = 1.0e-8

            solution = solve(
                prob, Rodas5P(); saveat = [0.0, 0.1], tolerances...
            )
            @test SciMLBase.successful_retcode(solution)
            @test solution.t[1] == 0.0
            @test solution.u[1] ≈ expected_state atol = 1.0e-8

            state_jacobian = ForwardDiff.jacobian(p) do dual_p
                dual_prob = make_override_dual_problem(
                    dual_p, inplace, vector_init_params
                )
                dual_solution = solve(
                    dual_prob, Rodas5P(); saveat = [0.0, 0.1], tolerances...
                )
                @test SciMLBase.successful_retcode(dual_solution)
                collect(dual_solution.u[1])
            end
            @test state_jacobian ≈ expected_jacobian atol = 1.0e-8

            residual = override_dual_init_residual(
                [solution.u[1][2]], (p[1], solution.u[1][1])
            )
            @test residual ≈ [0.0] atol = 1.0e-8
            residual_jacobian = ForwardDiff.jacobian(p) do dual_p
                dual_prob = make_override_dual_problem(
                    dual_p, inplace, vector_init_params
                )
                dual_solution = solve(
                    dual_prob, Rodas5P(); saveat = [0.0], tolerances...
                )
                @test SciMLBase.successful_retcode(dual_solution)
                c, _ = dual_p
                x, z = dual_solution.u[1]
                override_dual_init_residual([z], (c, x))
            end
            @test residual_jacobian ≈ zeros(1, 2) atol = 1.0e-8
        end
    end
end

function dual_rejecting_init_residual(u, p)
    eltype(u) === Float64 || throw(ArgumentError("initialization residual received Dual values"))
    return [u[1]^2 - 1]
end

@testset "default initialization uses finite differences when requested" begin
    initprob = NonlinearProblem(dual_rejecting_init_residual, [0.75])
    nlsolve = default_nlsolve(nothing, Val(false), initprob.u0, initprob, false)
    solution = solve(initprob, nlsolve)
    @test SciMLBase.successful_retcode(solution)
    @test solution.u[1] ≈ 1.0 atol = 1.0e-8
end
