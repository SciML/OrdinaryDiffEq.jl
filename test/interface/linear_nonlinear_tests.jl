using OrdinaryDiffEq, Test, Random, LinearAlgebra, LinearSolve, ADTypes
Random.seed!(123)

A = 0.01 * rand(3, 3)
rn = (du, u, p, t) -> begin
    mul!(du, A, u)
end
u0 = rand(3)
prob = ODEProblem(rn, u0, (0, 50.0))

function precsl(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pl = lu(convert(AbstractMatrix, W), check = false)
    else
        Pl = Plprev
    end
    return Pl, nothing
end

function precsr(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pr = lu(convert(AbstractMatrix, W), check = false)
    else
        Pr = Prprev
    end
    return nothing, Pr
end

function precslr(W, du, u, p, t, newW, Plprev, Prprev, solverdata)
    if newW === nothing || newW
        Pr = lu(convert(AbstractMatrix, W), check = false)
    else
        Pr = Prprev
    end
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
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precsl, smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precsr, smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precslr, smooth_est = false, concrete_jac = true
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
        linsolve = KrylovJL_GMRES(),
        precs = precslr, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    Rodas4(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precslr, concrete_jac = true
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
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precsl, smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precsr, smooth_est = false, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    TRBDF2(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precslr, smooth_est = false, concrete_jac = true
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
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precslr, concrete_jac = true
    )
);
@test length(sol.t) < 20
sol = @test_nowarn solve(
    prob,
    Rodas4(
        autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(),
        precs = precslr, concrete_jac = true
    )
);
@test length(sol.t) < 20
