using OrdinaryDiffEq, StochasticDiffEq, Test, Sundials

f_oop(u, p, t) = 2u
u0 = 0.5
tspan = (0.0, 1.0)
prob = ODEProblem(f_oop, u0, tspan)
sol = solve(prob, Tsit5())

function f_iip(du, u, p, t)
    return du .= 2.0 * u
end
prob = ODEProblem(f_iip, u0, tspan)
@test_throws SciMLBase.IncompatibleInitialConditionError solve(prob, Tsit5())

prob = ODEProblem{false}(f_oop, u0, tspan)
sol = solve(prob, Tsit5())
sol = solve(prob, nothing, alg = Tsit5())
sol = init(prob, nothing, alg = Tsit5())

prob = ODEProblem{false}(f_oop, 1.0 + im, tspan)
@test_throws SciMLBase.ComplexSupportError solve(prob, CVODE_Adams())

@test_throws SciMLBase.ProblemSolverPairingError solve(prob, DFBDF())
@test_throws SciMLBase.NonSolverError solve(prob, 5.0)

prob = ODEProblem{false}(f_oop, u0, (nothing, nothing))
@test_throws SciMLBase.NoTspanError solve(prob, Tsit5())

prob = ODEProblem{false}(f_oop, u0, (NaN, 1.0))
@test_throws SciMLBase.NaNTspanError solve(prob, Tsit5())

prob = ODEProblem{false}(f_oop, u0, (1.0, NaN))
@test_throws SciMLBase.NaNTspanError solve(prob, Tsit5())

prob = ODEProblem{false}(f_oop, Any[1.0, 1.0f0], tspan)
@test_throws SciMLBase.NonConcreteEltypeError solve(prob, Tsit5())

prob = ODEProblem{false}(f_oop, (1.0, 1.0f0), tspan)
@test_throws SciMLBase.TupleStateError solve(prob, Tsit5())

prob = ODEProblem{false}(f_oop, u0, (0.0 + im, 1.0))
@test_throws SciMLBase.ComplexTspanError solve(prob, Tsit5())

for u0 in ([0.0, 0.0], nothing)
    local fmm = ODEFunction(f_oop, mass_matrix = zeros(3, 3))
    local prob = ODEProblem(fmm, u0, (0.0, 1.0))
    @test_throws SciMLBase.IncompatibleMassMatrixError solve(prob, Tsit5())
end

# Allow empty mass matrix for empty u0
fmm = ODEFunction((du, u, t) -> nothing, mass_matrix = zeros(0, 0))
prob = ODEProblem(fmm, nothing, (0.0, 1.0))
sol = solve(prob, Tsit5())
@test isa(sol, SciMLBase.ODESolution)

f_sde(du, u, p, t) = du .= 1.01u
function g(du, u, p, t)
    du[1, 1] = 0.3u[1]
    du[1, 2] = 0.6u[1]
    du[1, 3] = 0.9u[1]
    du[1, 4] = 0.12u[1]
    du[2, 1] = 1.2u[2]
    du[2, 2] = 0.2u[2]
    du[2, 3] = 0.3u[2]
    return du[2, 4] = 1.8u[2]
end

prob = SDEProblem(
    f_sde,
    g,
    randn(ComplexF64, 2),
    (0.0, 1.0),
    noise_rate_prototype = complex(zeros(2, 4)),
    noise = StochasticDiffEq.RealWienerProcess(0.0, zeros(3))
)
@test_throws SciMLBase.NoiseSizeIncompatibilityError solve(prob, LambaEM())

function g!(du, u, p, t)
    du[1] .= u[1] + ones(3, 3)
    return du[2] .= ones(3, 3)
end
u0 = [zeros(3, 3), zeros(3, 3)]
prob = ODEProblem(g!, u0, (0, 1.0))
@test_throws SciMLBase.NonNumberEltypeError solve(prob, Tsit5())
