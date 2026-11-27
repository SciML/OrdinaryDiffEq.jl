using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqSDIRK, OrdinaryDiffEqNonlinearSolve,
    OrdinaryDiffEqBDF, StaticArrays, LinearAlgebra, SparseArrays, Test, ADTypes, ForwardDiff,
    SciMLBase
using OrdinaryDiffEqCore: default_nlsolve
using SymbolicIndexingInterface: parameter_values
import OrdinaryDiffEqCore, SymbolicIndexingInterface

if VERSION >= v"1.11"
    @test Base.ispublic(OrdinaryDiffEqCore, :default_nlsolve)
    @test Base.ispublic(SymbolicIndexingInterface, :parameter_values)
end

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
    integ = init(prob, Rodas5P())
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
        c2, _ = parameter_values(valp)
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

log_init_residual(u, p) = [log(u[1])]
log_init_jacobian(u, p) = reshape([inv(u[1])], 1, 1)
log_dae_rhs(u, p, t) = [-u[1], log(u[2])]

@testset "OverrideInit default solver on a logarithmic constraint" begin
    for analytic_jacobian in (false, true),
            alg in (Rodas5P(), Rodas5P(autodiff = AutoFiniteDiff()))

        nlf = analytic_jacobian ?
            NonlinearFunction(log_init_residual; jac = log_init_jacobian) :
            log_init_residual
        initprob = NonlinearProblem(nlf, [3.0])
        initialization_data = SciMLBase.OverrideInitData(
            initprob, nothing, nlsol -> [1.0, nlsol.u[1]], nothing
        )
        f = ODEFunction(
            log_dae_rhs; mass_matrix = Diagonal([1.0, 0.0]), initialization_data
        )
        prob = ODEProblem(f, [1.0, 3.0], (0.0, 0.1))
        sol = solve(prob, alg; abstol = 1.0e-10, reltol = 1.0e-10)
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[1] ≈ [1.0, 1.0] atol = 1.0e-8
        @test abs(log(sol.u[1][2])) < 1.0e-8
    end
end

function padded_init_residual!(resid, z, p)
    c, x = p
    resid[1] = z[1]^3 + z[1] - c * x
    for i in 2:length(z)
        resid[i] = z[i]
    end
    return nothing
end

function padded_dae_rhs!(du, u, p, t)
    c, k = p
    du[1] = -k * u[1] + u[2]
    du[2] = u[2]^3 + u[2] - c * u[1]
    for i in 3:length(u)
        du[i] = u[i]
    end
    return nothing
end

padded_init_guess(c, guess, n) = vcat(one(c) * guess, fill(zero(c), n - 1))
padded_init_params(c, named) = named ? (; c, x = c / 2) : (c, c / 2)

function padded_init_problem(c, n, guess, named)
    initprob = NonlinearProblem(
        padded_init_residual!, padded_init_guess(c, guess, n), padded_init_params(c, named)
    )
    update_initializeprob! = function (ip, valp)
        c2 = parameter_values(valp)[1]
        return SciMLBase.remake(
            ip; u0 = padded_init_guess(c2, guess, n), p = padded_init_params(c2, named)
        )
    end
    initialization_data = SciMLBase.OverrideInitData(
        initprob, update_initializeprob!, nlsol -> [nlsol.prob.p[2]; nlsol.u], nothing;
        is_update_oop = Val(true)
    )
    f = ODEFunction(
        padded_dae_rhs!; mass_matrix = Diagonal([1.0; zeros(n)]), initialization_data
    )
    return ODEProblem(f, [one(c); padded_init_guess(c, guess, n)], (0.0, 0.1), [c, 0.5])
end

# z^3 + z = c^2 / 2 at c = 2 has root z = 1 with dz/dc = c / (3z^2 + 1) = 0.5.
@testset "OverrideInit $(named ? "NamedTuple" : "Tuple") parameters, n = $n, guess = $guess" for (n, guess) in ((26, 0.5), (1, 1.0)),
        named in (false, true)

    tolerances = (; abstol = 1.0e-12, reltol = 1.0e-12)
    init_state(c) = init(padded_init_problem(c, n, guess, named), Rodas5P(); tolerances...).u
    saved_state(c) = solve(
        padded_init_problem(c, n, guess, named), Rodas5P(); saveat = [0.0], tolerances...
    ).u[1]
    @test init_state(2.0)[1:2] ≈ [1.0, 1.0] atol = 1.0e-8
    @test ForwardDiff.derivative(c -> init_state(c)[2], 2.0) ≈ 0.5 atol = 1.0e-8
    @test ForwardDiff.derivative(c -> saved_state(c)[2], 2.0) ≈ 0.5 atol = 1.0e-8
    residual_derivative = ForwardDiff.derivative(2.0) do c
        x, z = saved_state(c)
        z^3 + z - c * x
    end
    @test residual_derivative ≈ 0.0 atol = 1.0e-8
end

const mixed_root_slope = 2.0^-30
mixed_init_residual(u, p) = [u[1] - 1.0 - mixed_root_slope * p.c]
mixed_dae_rhs(u, p, t) = [-u[1], u[2] - 1.0 - mixed_root_slope * p.c]

function mixed_precision_problem(c, guess)
    initprob = NonlinearProblem(mixed_init_residual, [guess], (; c))
    initialization_data = SciMLBase.OverrideInitData(
        initprob, nothing, nlsol -> [1.0, nlsol.u[1]], nothing
    )
    f = ODEFunction(
        mixed_dae_rhs; mass_matrix = Diagonal([1.0, 0.0]), initialization_data
    )
    return ODEProblem(f, [1.0, guess], (0.0, 0.1), (; c))
end

# Float64 state, Float32 parameter seed: z = 1 + 2^-30 c has dz/dc = 2^-30, and the root
# 1 + 2^-29 is not representable in Float32.
@testset "OverrideInit keeps a Float64 root under a Float32 parameter seed" begin
    tolerances = (; abstol = 1.0e-12, reltol = 1.0e-12)
    root = 1.0 + mixed_root_slope * 2
    Tag = typeof(ForwardDiff.Tag(mixed_init_residual, Float32))
    seed = ForwardDiff.Dual{Tag}(2.0f0, 1.0f0)
    for guess in (0.5, root)
        initprob = NonlinearProblem(mixed_init_residual, [guess], (; c = seed))
        nlsol = solve(initprob, default_nlsolve(nothing, Val(false), [guess], initprob, true))
        @test SciMLBase.successful_retcode(nlsol)
        @test eltype(nlsol.u) === ForwardDiff.Dual{Tag, Float64, 1}
        @test ForwardDiff.value(nlsol.u[1]) ≈ root atol = tolerances.abstol rtol = 0
        @test ForwardDiff.partials(nlsol.u[1], 1) ≈ mixed_root_slope rtol = 1.0e-12

        integrator = init(mixed_precision_problem(seed, guess), Rodas5P(); tolerances...)
        sol = solve(mixed_precision_problem(seed, guess), Rodas5P(); tolerances...)
        @test SciMLBase.successful_retcode(sol)
        for z in (integrator.u[2], sol.u[1][2])
            @test ForwardDiff.value(z) isa Float64
            @test ForwardDiff.value(z) ≈ root atol = tolerances.abstol rtol = 0
            @test ForwardDiff.partials(z, 1) ≈ mixed_root_slope rtol = 1.0e-12
        end
    end
end

two_param_init_residual(z, p) = [z[1]^3 + z[1] - p[1] * p[2]]
two_param_dae_rhs(u, p, t) = [-u[1], u[2]^3 + u[2] - p[1] * p[2]]

function two_param_problem(a, b, named)
    params = named ? (; a, b) : (a, b)
    initprob = NonlinearProblem(two_param_init_residual, [0.5], params)
    initialization_data = SciMLBase.OverrideInitData(
        initprob, nothing, nlsol -> [1.0, nlsol.u[1]], nothing
    )
    f = ODEFunction(
        two_param_dae_rhs; mass_matrix = Diagonal([1.0, 0.0]), initialization_data
    )
    return ODEProblem(f, [1.0, 0.5], (0.0, 0.1), params)
end

# z^3 + z = a b at (a, b) = (2, 1) has root z = 1 and ∇z = (b, a) / (3z^2 + 1) = (1/4, 1/2).
@testset "OverrideInit partials for two seeded $(named ? "NamedTuple" : "Tuple") parameters" for named in (false, true)
    tolerances = (; abstol = 1.0e-12, reltol = 1.0e-12)
    gradient = ForwardDiff.gradient([2.0, 1.0]) do ab
        sol = solve(two_param_problem(ab..., named), Rodas5P(); saveat = [0.0], tolerances...)
        @test SciMLBase.successful_retcode(sol)
        sol.u[1][2]
    end
    @test gradient ≈ [0.25, 0.5] atol = 1.0e-10
end

const override_root_seed = ForwardDiff.Dual{typeof(ForwardDiff.Tag(vcat, Float64))}(2.0, 1.0)

# The direct nonlinear solve, the `init` state and the first saved state for an in-place DAE
# whose algebraic states are initialized by `nlf` with parameters `(; c = override_root_seed)`.
function override_root_solutions(
        nlf, z0, dae_residual; autodiff = true, iip_selector = Val(false)
    )
    tolerances = (; abstol = 1.0e-12, reltol = 1.0e-12)
    c = override_root_seed
    nlprob = NonlinearProblem(nlf, z0, (; c))
    nlsolve = default_nlsolve(nothing, iip_selector, z0, nlprob, autodiff)
    nlsol = solve(nlprob, nlsolve; tolerances...)
    initialization_data = SciMLBase.OverrideInitData(
        nlprob, nothing, sol -> vcat(1.0, sol.u), nothing
    )
    rhs! = function (du, u, p, t)
        du[1] = -u[1]
        du[2:end] .= dae_residual(u[2:end], p)
        return nothing
    end
    f = ODEFunction{true}(
        rhs!; mass_matrix = Diagonal([1.0; zeros(length(z0))]), initialization_data
    )
    prob = ODEProblem(f, Vector(vcat(one(c), z0)), (0.0, 0.1), (; c))
    alg = autodiff ? Rodas5P() : Rodas5P(autodiff = AutoFiniteDiff())
    integrator = init(prob, alg; tolerances...)
    sol = solve(prob, alg; saveat = [0.0], tolerances...)
    @test SciMLBase.successful_retcode(nlsol)
    @test SciMLBase.successful_retcode(sol)
    return [vcat(nlsol.u), integrator.u[2:end], sol.u[1][2:end]]
end

function test_override_root(states, root, tangent)
    for z in states
        @test ForwardDiff.value.(z) ≈ root atol = 1.0e-8
        @test ForwardDiff.partials.(z, 1) ≈ tangent atol = 1.0e-8
    end
    return nothing
end

# (z - c)^3 = 0 has a singular Jacobian at its root z = c, so a consistent dual warm start
# z = c keeps its tangent 1. Forward differences estimate that Jacobian as about 9e-16.
@testset "OverrideInit repeated $label root keeps a consistent dual warm start, autodiff = $autodiff" for (label, residual, z0) in (
            ("scalar", (z, p) -> (z - p.c)^3, override_root_seed),
            ("vector", (z, p) -> [(z[1] - p.c)^3], [override_root_seed]),
        ),
        autodiff in (true, false)

    states = override_root_solutions(residual, z0, (z, p) -> z .- p.c; autodiff)
    test_override_root(states, [2.0], [1.0])
end

# 2^-30 (z - c) = 0 has the regular root z = c with tangent 1. Its finite-difference
# Jacobian is resolved, so the partials come from implicit differentiation. A residual that
# sees only the value of z has a ForwardDiff Jacobian of zero, so it is checked with
# finite differences only.
small_root_residual(z, p) = [2.0^-30 * (z[1] - p.c)]
small_root_value_residual(z, p) = [2.0^-30 * (ForwardDiff.value(z[1]) - p.c)]
@testset "OverrideInit small regular Jacobian, $label residual, autodiff = $autodiff" for (label, residual, autodiff) in (
        ("dual", small_root_residual, true),
        ("dual", small_root_residual, false),
        ("value-only", small_root_value_residual, false),
    )

    nlprob = NonlinearProblem(residual, [0.5], (; c = override_root_seed))
    nlsolve = default_nlsolve(nothing, Val(false), [0.5], nlprob, autodiff)
    @test solve(nlprob, nlsolve; abstol = 1.0e-12, reltol = 1.0e-12).alg === nlsolve
    states = override_root_solutions(
        residual, [0.5], (z, p) -> 2.0^-30 .* (z .- p.c); autodiff
    )
    test_override_root(states, [2.0], [1.0])
end

# 2z₁ + z₂ = c, z₁ + 3z₂ = 0 has root (3c, -c) / 5 and tangent (0.6, -0.2).
@testset "OverrideInit with a sparse analytic Jacobian" begin
    residual(z, p) = [2z[1] + z[2] - p.c, z[1] + 3z[2]]
    nlf = NonlinearFunction(residual; jac = (z, p) -> sparse([2.0 1.0; 1.0 3.0]))
    states = override_root_solutions(nlf, [0.0, 0.0], residual)
    test_override_root(states, [1.2, -0.4], [0.6, -0.2])
end

@testset "OverrideInit of an in-place DAE from a static nonlinear warm start" begin
    nlf = NonlinearFunction(
        (z, p) -> SVector(z[1] - p.c); jac = (z, p) -> SMatrix{1, 1}(1.0)
    )
    states = override_root_solutions(
        nlf, SVector(override_root_seed), (z, p) -> z .- p.c; iip_selector = Val(true)
    )
    @test states[1] isa SVector
    test_override_root(states, [2.0], [1.0])
end

# The state enters the initialization residual only through its value, so the Jacobian
# must come from the finite-difference backend, giving z = c / 2 with tangent 0.5.
@testset "OverrideInit with AutoFiniteDiff, analytic Jacobian = $analytic" for analytic in (false, true)
    jac = analytic ? ((z, p) -> fill(2.0, 1, 1)) : nothing
    nlf = NonlinearFunction((z, p) -> [2 * ForwardDiff.value(z[1]) - p.c]; jac)
    states = override_root_solutions(
        nlf, [0.5], (z, p) -> 2 .* z .- p.c; autodiff = false
    )
    test_override_root(states, [1.0], [0.5])
    for z in states
        @test ForwardDiff.partials(2 * z[1] - override_root_seed, 1) ≈ 0.0 atol = 1.0e-8
    end
end

# records the keywords Brown passes to the nonlinear solve
struct SpyNLSolve <: SciMLBase.AbstractNonlinearAlgorithm
    kwargs::Base.RefValue{Any}
end
SpyNLSolve() = SpyNLSolve(Ref{Any}((;)))
function SciMLBase.__solve(prob::NonlinearProblem, alg::SpyNLSolve, args...; kwargs...)
    alg.kwargs[] = NamedTuple(kwargs)
    return SciMLBase.build_solution(prob, alg, prob.u0, prob.u0; retcode = ReturnCode.Success)
end

@testset "BrownFullBasicInit tolerances" begin
    # y₁' = -y₁, 0 = y₂³ + y₂ - y₁
    mm_oop(u, p, t) = [-u[1], u[2]^3 + u[2] - u[1]]
    mm_iip(du, u, p, t) = (du .= mm_oop(u, p, t); nothing)
    dae_oop(du, u, p, t) = [du[1] + u[1], u[2]^3 + u[2] - u[1]]
    dae_iip(r, du, u, p, t) = (r .= dae_oop(du, u, p, t); nothing)
    y₂ = 0.6823278038280193 # root of y₂³ + y₂ = 1

    M = Diagonal([1.0, 0.0])
    dv = [true, false]
    make_probs(u0; kw...) = [
        (ODEProblem(ODEFunction(mm_oop; mass_matrix = M), u0, (0.0, 1.0); kw...), FBDF()),
        (ODEProblem(ODEFunction(mm_iip; mass_matrix = M), u0, (0.0, 1.0); kw...), FBDF()),
        (DAEProblem(dae_oop, [-1.0, 0.0], u0, (0.0, 1.0); differential_vars = dv, kw...), DFBDF()),
        (DAEProblem(dae_iip, [-1.0, 0.0], u0, (0.0, 1.0); differential_vars = dv, kw...), DFBDF()),
    ]

    # per-state tolerances are reduced to scalars for the nonlinear solve
    for (prob, alg) in make_probs([1.0, 2.0])
        for initializealg in (BrownFullBasicInit(), BrownFullBasicInit(abstol = nothing))
            integ = init(prob, alg; initializealg, abstol = [1.0e-8, 1.0e-8], reltol = [1.0e-8, 1.0e-8])
            @test integ.u ≈ [1.0, y₂]
        end
    end

    # the residual is below the solver's abstol but above Brown's default of 1e-10,
    # so only `abstol = nothing` leaves u0 untouched
    u0 = [1.0, y₂ + 1.0e-5]
    for (prob, alg) in make_probs(u0), abstol in (1.0e-3, [1.0e-3, 1.0e-3])
        integ = init(prob, alg; initializealg = BrownFullBasicInit(abstol = nothing), abstol)
        @test integ.u == u0
        integ = init(prob, alg; initializealg = BrownFullBasicInit(), abstol)
        @test integ.u ≈ [1.0, y₂]
        @test integ.u != u0
    end

    # the tolerances reach the nonlinear solve, with vectors reduced to their tightest entry
    for (prob, alg) in make_probs([1.0, 2.0])
        spy = SpyNLSolve()
        init(
            prob, alg; initializealg = BrownFullBasicInit(; abstol = 1.0e-6, nlsolve = spy),
            abstol = 1.0e-3, reltol = [1.0e-4, 1.0e-5]
        )
        @test get(spy.kwargs[], :abstol, nothing) == 1.0e-6
        @test get(spy.kwargs[], :reltol, nothing) == 1.0e-5

        spy = SpyNLSolve()
        init(
            prob, alg; initializealg = BrownFullBasicInit(; abstol = nothing, nlsolve = spy),
            abstol = [1.0e-2, 1.0e-3]
        )
        @test get(spy.kwargs[], :abstol, nothing) == 1.0e-3
    end

    # `abstol = nothing` also works when the algorithm is set on the problem
    spy = SpyNLSolve()
    initializealg = BrownFullBasicInit(; abstol = nothing, nlsolve = spy)
    for (prob, alg) in make_probs([1.0, 2.0]; initializealg)
        spy.kwargs[] = (;)
        init(prob, alg; abstol = 1.0e-3)
        @test get(spy.kwargs[], :abstol, nothing) == 1.0e-3
    end
end

@testset "BrownFullBasicInit with out-of-place static arrays" begin
    f_static(u, p, t) = SA[-u[1] + u[2], u[1] + u[2] - p]
    fn = ODEFunction(f_static, mass_matrix = [1.0 0.0; 0.0 0.0])
    prob = ODEProblem(fn, SA[1.0, 0.5], (0.0, 1.0), 1.0)
    sol = solve(prob, Rodas5P(); initializealg = BrownFullBasicInit())
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[1] isa SVector{2, Float64}
    @test sol.u[1] ≈ SA[1.0, 0.0] atol = 1.0e-10

    # Reinitialization after a callback changes the algebraic variable only.
    cb = DiscreteCallback(
        (u, t, integ) -> t == 0.5, integ -> (integ.p = 2.0);
        initializealg = BrownFullBasicInit()
    )
    sol = solve(
        prob, Rodas5P(); initializealg = BrownFullBasicInit(), callback = cb,
        tstops = [0.5], abstol = 1.0e-10, reltol = 1.0e-10
    )
    @test SciMLBase.successful_retcode(sol)
    i = findall(==(0.5), sol.t)
    @test length(i) == 2
    @test sol.u[i[2]][1] == sol.u[i[1]][1]
    @test sol.u[i[2]][2] ≈ 2.0 - sol.u[i[1]][1] atol = 1.0e-10
    @test sol.u[i[1]][2] ≈ 1.0 - sol.u[i[1]][1] atol = 1.0e-8
end

@testset "BrownFullBasicInit out-of-place mutable u0 aliasing" begin
    f_vec(u, p, t) = [-u[1] + u[2], u[1] + u[2] - 1.0]
    fn = ODEFunction(f_vec, mass_matrix = [1.0 0.0; 0.0 0.0])

    prob = ODEProblem(fn, [1.0, 0.5], (0.0, 1.0))
    integ = init(
        prob, Rodas5P(); initializealg = BrownFullBasicInit(),
        alias = SciMLBase.ODEAliasSpecifier(alias_u0 = true)
    )
    @test integ.u === prob.u0
    @test prob.u0 ≈ [1.0, 0.0] atol = 1.0e-10

    prob = ODEProblem(fn, [1.0, 0.5], (0.0, 1.0))
    integ = init(prob, Rodas5P(); initializealg = BrownFullBasicInit())
    @test integ.u !== prob.u0
    @test prob.u0 == [1.0, 0.5]
    @test integ.u ≈ [1.0, 0.0] atol = 1.0e-10
end
