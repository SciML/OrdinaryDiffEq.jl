using OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF, OrdinaryDiffEqRosenbrock, Test, Random,
    LinearAlgebra, LinearSolve, ADTypes, SciMLBase
using OrdinaryDiffEqNonlinearSolve: NonlinearSolveAlg
using NonlinearSolve: NewtonRaphson
using SimpleNonlinearSolve: SimpleNewtonRaphson
Random.seed!(123)

A = 0.01 * rand(3, 3)
rn = (du, u, p, t) -> begin
    mul!(du, A, u)
end
u0 = rand(3)
prob = ODEProblem(rn, u0, (0, 50.0))

# Preconditioners for the new interface: precs(A, p) -> (Pl, Pr)
# where p = (du, u, params, t) is passed through the LinearProblem
function precsl(A, p)
    Pl = lu(convert(AbstractMatrix, A), check = false)
    return Pl, I
end

function precsr(A, p)
    Pr = lu(convert(AbstractMatrix, A), check = false)
    return I, Pr
end

function precslr(A, p)
    Pr = lu(convert(AbstractMatrix, A), check = false)
    return Pr, Pr
end

sol = @test_nowarn solve(prob, TRBDF2(autodiff = AutoFiniteDiff()));
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES())
);
@test length(sol.t) < 20
solref = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        smooth_est = false
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precsl),
        smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precsr),
        smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precslr),
        smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    QNDF(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        concrete_jac = true
    )
);
@test length(sol.t) < 25
sol = @test_nowarn solve(
    prob,
    Rosenbrock23(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precslr),
        concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    Rodas4(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precslr),
        concrete_jac = true
    )
);
@test length(sol.t) < 20

sol = @test_nowarn solve(prob, TRBDF2(autodiff = AutoFiniteDiff()));
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob, TRBDF2(autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES())
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        smooth_est = false
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precsl),
        smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precsr),
        smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precslr),
        smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    QNDF(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        concrete_jac = true
    )
);
@test length(sol.t) < 25
sol = @test_nowarn solve(
    prob,
    Rosenbrock23(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precslr),
        concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    Rodas4(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(precs = precslr),
        concrete_jac = true
    )
);
@test length(sol.t) < 20

sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(),
        nlsolve = NonlinearSolveAlg(NewtonRaphson())
    )
);
@test length(sol.t) < 20
let integ = init(
        prob,
        TRBDF2(
            autodiff = AutoFiniteDiff(),
            nlsolve = NonlinearSolveAlg(NewtonRaphson())
        )
    )
    @test integ.cache.nlsolver.cache.weight !== nothing
    step!(integ)
    @test !iszero(integ.cache.nlsolver.cache.weight)
end

sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(),
        linsolve = KrylovJL_GMRES(),
        nlsolve = NonlinearSolveAlg(NewtonRaphson(autodiff = AutoFiniteDiff())),
        concrete_jac = true
    )
);
@test length(sol.t) < 20
let integ = init(
        prob,
        TRBDF2(
            autodiff = AutoFiniteDiff(),
            linsolve = KrylovJL_GMRES(),
            nlsolve = NonlinearSolveAlg(NewtonRaphson(autodiff = AutoFiniteDiff())),
            concrete_jac = true
        )
    )
    @test integ.cache.nlsolver.cache.weight !== nothing
    step!(integ)
    @test !iszero(integ.cache.nlsolver.cache.weight)
end

using StaticArrays
let
    vdp_static(u, p, t) = SVector(u[2], p[1] * ((1 - u[1]^2) * u[2] - u[1]))
    prob = ODEProblem(
        vdp_static, SVector(2.0, 0.0), (0.0, 1.0), SVector(1.0e3)
    )
    integ = init(
        prob,
        TRBDF2(
            nlsolve = NonlinearSolveAlg(SimpleNewtonRaphson(; autodiff = AutoForwardDiff())),
            concrete_jac = true
        ); reltol = 1.0e-8, abstol = 1.0e-10
    )
    @test integ.cache.nlsolver.cache.W isa Base.RefValue
    sol = solve!(integ)
    @test SciMLBase.successful_retcode(sol.retcode)
    solref = solve(prob, TRBDF2(); reltol = 1.0e-8, abstol = 1.0e-10)
    @test maximum(abs.(sol.u[end] .- solref.u[end])) < 1.0e-3
end
