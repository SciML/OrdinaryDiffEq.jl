using OrdinaryDiffEq, Base.Test, DiffEqOperators
N = 100
dx = 1.0; dt=0.01
srand(0); u0 = rand(N)
reltol = 1e-4
# L = DerivativeOperator{Float64}(2,2,dx,N,:Dirichlet0,:Dirichlet0) # error for caching version at the moment
dd = -2.0 * ones(N); du = 1.0 * ones(N-1)
A = diagm(du,-1) + diagm(dd,0) + diagm(du,1); L = DiffEqArrayOperator(A)
krylov_f2 = (u,p,t) -> -0.1*u
krylov_f2! = (du,u,p,t) -> du .= -0.1*u
prob = SplitODEProblem(L,krylov_f2,u0,(0.0,1.0))
prob_inplace = SplitODEProblem(L,krylov_f2!,u0,(0.0,1.0))

# Ad-hoc fix for SplitFunction miscalssified as having analytic solutions
DiffEqBase.has_analytic(::typeof(prob.f)) = false
DiffEqBase.has_analytic(::typeof(prob_inplace.f)) = false

sol = solve(prob, LawsonEuler(); dt=dt)
sol_krylov = solve(prob, LawsonEuler(krylov=true); dt=dt, reltol=reltol)
@test isapprox(sol.u,sol_krylov.u; rtol=reltol)

sol_ip = solve(prob_inplace, LawsonEuler(); dt=0.01)
sol_ip_krylov = solve(prob_inplace, LawsonEuler(krylov=true); dt=dt, reltol=reltol)
@test isapprox(sol_ip.u,sol_ip_krylov.u; rtol=reltol)
