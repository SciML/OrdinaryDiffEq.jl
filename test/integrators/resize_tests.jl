using OrdinaryDiffEq, Test
f(du, u, p, t) = du .= u
prob = ODEProblem(f, [1.0], (0.0, 1.0))

i = init(prob, Tsit5())
resize!(i, 5)
@test length(i.cache.u) == 5
@test length(i.cache.uprev) == 5
@test length(i.cache.k1) == 5
@test length(i.cache.k2) == 5
@test length(i.cache.k3) == 5
@test length(i.cache.k4) == 5
@test length(i.cache.k5) == 5
@test length(i.cache.k6) == 5
@test length(i.cache.k7) == 5
solve!(i)

i = init(prob, ImplicitEuler())
resize!(i, 5)
@test length(i.cache.atmp) == 5
@test length(i.cache.uprev) == 5
# nlsolver fields
@test length(i.cache.nlsolver.z) == 5
@test length(i.cache.nlsolver.dz) == 5
@test length(i.cache.nlsolver.weight) == 5
@test length(i.cache.nlsolver.ztmp) == 5
@test length(i.cache.nlsolver.tmp) == 5
@test length(i.cache.nlsolver.k) == 5
@test length(i.cache.nlsolver.du1) == 5
# ForwardDiff
@test length(i.cache.nlsolver.jac_config.duals[1]) == 5
@test length(i.cache.nlsolver.jac_config.duals[2]) == 5
@test size(i.cache.nlsolver.cache.W) == (5,5)
@test size(i.cache.nlsolver.cache.J) == (5,5)
solve!(i)

i = init(prob, ImplicitEuler(;autodiff=false))
resize!(i, 5)
@test length(i.cache.atmp) == 5
@test length(i.cache.uprev) == 5
# nlsolver fields
@test length(i.cache.nlsolver.z) == 5
@test length(i.cache.nlsolver.dz) == 5
@test length(i.cache.nlsolver.weight) == 5
@test length(i.cache.nlsolver.ztmp) == 5
@test length(i.cache.nlsolver.tmp) == 5
@test length(i.cache.nlsolver.k) == 5
@test length(i.cache.nlsolver.du1) == 5
# DiffEqDiffTools
@test length(i.cache.nlsolver.jac_config.x1) == 5
@test length(i.cache.nlsolver.jac_config.fx) == 5
@test length(i.cache.nlsolver.jac_config.fx1) == 5
@test size(i.cache.nlsolver.cache.W) == (5,5)
@test size(i.cache.nlsolver.cache.J) == (5,5)
solve!(i)

i = init(prob, Rosenbrock23())
resize!(i, 5)
@test length(i.cache.u) == 5
@test length(i.cache.uprev) == 5
@test length(i.cache.k₁) == 5
@test length(i.cache.k₂) == 5
@test length(i.cache.k₃) == 5
@test length(i.cache.du1) == 5
@test length(i.cache.du2) == 5
@test length(i.cache.f₁) == 5
@test length(i.cache.fsalfirst) == 5
@test length(i.cache.fsallast) == 5
@test length(i.cache.dT) == 5
@test length(i.cache.tmp) == 5
@test size(i.cache.J) == (5, 5)
@test size(i.cache.W) == (5, 5)
@test length(i.cache.linsolve_tmp) == 5
@test length(i.cache.jac_config.duals[1]) == 5
@test length(i.cache.jac_config.duals[2]) == 5
solve!(i)

i = init(prob, Rosenbrock23(;autodiff=false))
resize!(i, 5)
@test length(i.cache.u) == 5
@test length(i.cache.uprev) == 5
@test length(i.cache.k₁) == 5
@test length(i.cache.k₂) == 5
@test length(i.cache.k₃) == 5
@test length(i.cache.du1) == 5
@test length(i.cache.du2) == 5
@test length(i.cache.f₁) == 5
@test length(i.cache.fsalfirst) == 5
@test length(i.cache.fsallast) == 5
@test length(i.cache.dT) == 5
@test length(i.cache.tmp) == 5
@test size(i.cache.J) == (5, 5)
@test size(i.cache.W) == (5, 5)
@test length(i.cache.linsolve_tmp) == 5
@test length(i.cache.jac_config.x1) == 5
@test length(i.cache.jac_config.fx) == 5
@test length(i.cache.jac_config.fx1) == 5
solve!(i)

function f(du,u,p,t)
  du[1] = 2.0 * u[1] - 1.2 * u[1]*u[2]
  du[2] = -3 * u[2] + u[1]*u[2]
  for i in 3:length(u)
  	du[i] = 0.0
  end
end
function f_jac(J,u,p,t)
  J[1,1] = 2.0 - 1.2 * u[2]
  J[1,2] = -1.2 * u[1]
  J[2,1] = 1 * u[2]
  J[2,2] = -3 + u[1]
  for i in 3:length(u)
  	for j in 3:length(u)
  		if i == j
  			J[i, j] = 1.0
  		else
	  	    J[i, j] = 0.0
	  	end
  	end
  end
  nothing
end
ff = ODEFunction(f;jac=f_jac,jac_prototype=[1.0 1.0; 1.0 1.0])

cb = DiscreteCallback((u,t,integ) -> true, integ -> @views(integ.u[3:5]) .= 0)
prob = ODEProblem(ff, [1.0, 1.0], (0.0, 1.0))
i = init(prob, ImplicitEuler(), callback=cb)
resize!(i, 5)
sol = solve!(i)
