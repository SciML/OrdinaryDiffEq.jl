using OrdinaryDiffEq, Base.Test, DiffEqOperators
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
    gc()
    sol = solve(prob, Alg(); dt=dt, internalnorm=Base.norm)
    sol_krylov = solve(prob, Alg(krylov=true, m=10); dt=dt, reltol=reltol, internalnorm=Base.norm)
    @test isapprox(sol.u,sol_krylov.u; rtol=reltol)

    sol_ip = solve(prob_inplace, Alg(); dt=dt, internalnorm=Base.norm)
    sol_ip_krylov = solve(prob_inplace, Alg(krylov=true, m=10); dt=dt, reltol=reltol, internalnorm=Base.norm)
    @test isapprox(sol.u,sol_krylov.u; rtol=reltol)

    println(Alg) # prevent Travis hanging
end
