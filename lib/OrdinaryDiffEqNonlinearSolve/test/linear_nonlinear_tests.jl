using OrdinaryDiffEqSDIRK, OrdinaryDiffEqBDF, OrdinaryDiffEqRosenbrock, Test, Random,
    LinearAlgebra, LinearSolve, ADTypes
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
