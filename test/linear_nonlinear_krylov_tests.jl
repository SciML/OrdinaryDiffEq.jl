using OrdinaryDiffEq, Test, DiffEqOperators
@testset "Classical ExpRK" begin
    N = 20
    dx = 1.0; dt=0.1
    srand(0); u0 = rand(N)
    reltol = 1e-4
    L = DerivativeOperator{Float64}(2,2,dx,N,:Dirichlet0,:Dirichlet0)
    krylov_f2 = (u,p,t) -> -0.1*u
    krylov_f2! = (du,u,p,t) -> du .= -0.1*u
    prob = SplitODEProblem(L,krylov_f2,u0,(0.0,1.0))
    prob_inplace = SplitODEProblem(L,krylov_f2!,u0,(0.0,1.0))

    # Ad-hoc fix for SplitFunction miscalssified as having analytic solutions
    DiffEqBase.has_analytic(::typeof(prob.f)) = false
    DiffEqBase.has_analytic(::typeof(prob_inplace.f)) = false

    Algs = [LawsonEuler,NorsettEuler,ETDRK2,ETDRK3,ETDRK4,HochOst4]
    for Alg in Algs
        sol = solve(prob, Alg(); dt=dt, internalnorm=Base.norm)
        sol_krylov = solve(prob, Alg(krylov=true, m=10); dt=dt, reltol=reltol, internalnorm=Base.norm)
        @test isapprox(sol.u,sol_krylov.u; rtol=reltol)

        sol_ip = solve(prob_inplace, Alg(); dt=dt, internalnorm=Base.norm)
        sol_ip_krylov = solve(prob_inplace, Alg(krylov=true, m=10); dt=dt, reltol=reltol, internalnorm=Base.norm)
        @test isapprox(sol.u,sol_krylov.u; rtol=reltol)

        println(Alg) # prevent Travis hanging
    end
end

@testset "EPIRK" begin
    N = 20
    srand(0); u0 = normalize(randn(N))
    # For the moment, use dense Jacobian
    dd = -2 * ones(N); du = ones(N-1)
    A = diagm(du, -1) + diagm(dd) + diagm(du, 1)
    _f = (u,p,t) -> A*u - u.^3
    _f_ip = (du,u,p,t) -> (mul!(du, A, u); du .-= u.^3)
    _jac = (u,p,t) -> A - 3 * diagm(u.^2)
    _jac_ip = (J,u,p,t) -> begin
        copyto!(J, A)
        @inbounds for i = 1:N
            J[i, i] -= 3 * u[i]^2
        end
    end
    f = ODEFunction{false}(_f; jac=_jac)
    f_ip = ODEFunction{true}(_f_ip; jac=_jac_ip)
    prob = ODEProblem(f, u0, (0.0, 1.0))
    prob_ip = ODEProblem{true}(f_ip, u0, (0.0, 1.0))

    dt = 0.05; tol=1e-5
    Algs = [Exp4, EPIRK4s3A, EPIRK4s3B, EXPRB53s3, EPIRK5P1, EPIRK5P2]
    for Alg in Algs
        sol = solve(prob, Alg(); dt=dt, internalnorm=Base.norm, reltol=tol)
        sol_ref = solve(prob, Tsit5(); reltol=tol)
        @test isapprox(sol(1.0), sol_ref(1.0); rtol=tol)

        sol = solve(prob_ip, Alg(); dt=dt, internalnorm=Base.norm, reltol=tol)
        sol_ref = solve(prob_ip, Tsit5(); reltol=tol)
        @test isapprox(sol(1.0), sol_ref(1.0); rtol=tol)
        println(Alg) # prevent Travis hanging
    end

    gc()
    sol = solve(prob, EPIRK5s3(); dt=dt, internalnorm=Base.norm, reltol=tol)
    sol_ref = solve(prob, Tsit5(); reltol=tol)
    @test_broken isapprox(sol(1.0), sol_ref(1.0); rtol=tol)

    sol = solve(prob_ip, EPIRK5s3(); dt=dt, internalnorm=Base.norm, reltol=tol)
    sol_ref = solve(prob_ip, Tsit5(); reltol=tol)
    @test_broken isapprox(sol(1.0), sol_ref(1.0); rtol=tol)
    println(EPIRK5s3) # prevent Travis hanging
end
