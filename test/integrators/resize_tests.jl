using OrdinaryDiffEq, Test
f(du, u, p, t) = du .= 0
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
@test length(i.cache.z) == 5
@test length(i.cache.dz) == 5
@test length(i.cache.b) == 5
@test length(i.cache.atmp) == 5
@test length(i.cache.tmp) == 5
@test length(i.cache.k) == 5
@test length(i.cache.du1) == 5
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
@test length(i.cache.jac_config.duals[1]) == 5
@test length(i.cache.jac_config.duals[2]) == 5
@test size(i.cache.J) == (5,5)
@test size(i.cache.W) == (5,5)
@test size(i.cache.nlsolver.cache.W) == (5,5)
solve!(i)

i = init(prob, ImplicitEuler(;autodiff=false))
resize!(i, 5)
@test length(i.cache.z) == 5
@test length(i.cache.dz) == 5
@test length(i.cache.b) == 5
@test length(i.cache.atmp) == 5
@test length(i.cache.tmp) == 5
@test length(i.cache.k) == 5
@test length(i.cache.du1) == 5
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
@test size(i.cache.J) == (5,5)
@test size(i.cache.W) == (5,5)
@test size(i.cache.nlsolver.cache.W) == (5,5)
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

# i = init(prob, Rosenbrock23(;autodiff=false))
# resize!(i, 5)
# @test length(i.cache.u) == 5
# @test length(i.cache.uprev) == 5
# @test length(i.cache.k₁) == 5
# @test length(i.cache.k₂) == 5
# @test length(i.cache.k₃) == 5
# @test length(i.cache.du1) == 5
# @test length(i.cache.du2) == 5
# @test length(i.cache.f₁) == 5
# @test length(i.cache.fsalfirst) == 5
# @test length(i.cache.fsallast) == 5
# @test length(i.cache.dT) == 5
# @test length(i.cache.tmp) == 5
# @test size(i.cache.J) == (5, 5)
# @test size(i.cache.W) == (5, 5)
# @test length(i.cache.linsolve_tmp) == 5
# @test length(i.cache.jac_config.x1) == 5
# @test length(i.cache.jac_config.fx) == 5
# @test length(i.cache.jac_config.fx1) == 5

function jac(J,u,p,t)
	l = length(u)
	for i in 1:l
		for j in 1:l
			if i == j
				J[i][j] = 1.0
			else
				J[i][j] = 0.0
			end
		end
	end
end

# prob = ODEProblem(ODEFunction(f;jac=(J,u,p,t)->J=jac,jac_prototype =[1.0 1.0; 1.0 1.0]), [1.0, 1.0], (0.0, 10.0))
# i = init(prob, ImplicitEuler())
# resize!(i, 5)
# @show i.cache.nlsolver.cache.W.J
# solve!(i)