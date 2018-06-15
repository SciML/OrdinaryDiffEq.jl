using OrdinaryDiffEq, Base.Test, DiffEqOperators
N = 100
dx = 1.0; dt=0.01
srand(0); u0 = rand(N)
reltol = 1e-4
L1 = DerivativeOperator{Float64}(2,4,dx,N,:Dirichlet0,:Dirichlet0)
L2 = DerivativeOperator{Float64}(4,4,dx,N,:Dirichlet0,:Dirichlet0)
L = 1.01*L1 + 2.02*L2
krylov_f2 = (u,p,t) -> -0.1*u
krylov_f2! = (du,u,p,t) -> du .= -0.1*u
prob = SplitODEProblem(L,krylov_f2,u0,(0.0,1.0))
prob_inplace = SplitODEProblem(L,krylov_f2!,u0,(0.0,1.0))

# Ad-hoc fix for SplitFunction miscalssified as having analytic solutions
DiffEqBase.has_analytic(::typeof(prob.f)) = false
DiffEqBase.has_analytic(::typeof(prob_inplace.f)) = false

Algs = [LawsonEuler,NorsettEuler,ExpTrapezoid]
for Alg in Algs
    gc()
    sol = solve(prob, Alg(); dt=dt, internalnorm=Base.norm)
    sol_krylov = solve(prob, Alg(krylov=true); dt=dt, reltol=reltol, internalnorm=Base.norm)
    @test isapprox(sol.u,sol_krylov.u; rtol=reltol)

    sol_ip = solve(prob_inplace, Alg(); dt=dt, internalnorm=Base.norm)
    sol_ip_krylov = solve(prob_inplace, Alg(krylov=true); dt=dt, reltol=reltol, internalnorm=Base.norm)
    @test isapprox(sol_ip.u,sol_ip_krylov.u; rtol=reltol)
end
