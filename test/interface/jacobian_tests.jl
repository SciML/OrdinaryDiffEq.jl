using OrdinaryDiffEq, ForwardDiff, Test

function d_alembert(du,u,p,t)
  du[1] = p[1] - p[2]*u[1] + p[3]*t
end

function d_alembert_jac(J,u,p,t)
  J[1] = -p[2]
end

function d_alembert_analytic(u0,p,t::Number)
    a,b,c = p
    ebt = exp(b*t)
    @. exp(-b*t)*(-a*b + c + ebt*(a*b + c*(b*t - 1)) + b^2 * u0)/(b^2)
end

p = (1., 2., 3.)
u0 = [1.0]
tspan = (0.0,10.0)
prob = ODEProblem(ODEFunction(d_alembert,
                              jac = d_alembert_jac,
                              analytic=d_alembert_analytic),
                              u0,tspan,p)

sol = solve(prob,Tsit5(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Rosenbrock23(),abstol=1e-8,reltol=1e-8)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Rodas4(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Veldd4(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,Rodas5(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,TRBDF2(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 2e-6
sol = solve(prob,Trapezoid(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 2e-6
sol = solve(prob,KenCarp3(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 8e-4
sol = solve(prob,KenCarp4(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,KenCarp47(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7
sol = solve(prob,KenCarp58(),abstol=1e-10,reltol=1e-10)
@test sol.errors[:l2] < 1e-7

using ModelingToolkit
function lotka(du,u,p,t)
  x = u[1]
  y = u[2]
  du[1] = p[1]*x - p[2]*x*y
  du[2] = -p[3]*y + p[4]*x*y
end

prob = ODEProblem(lotka,[1.0,1.0],(0.0,1.0),[1.5,1.0,3.0,1.0])
de = ModelingToolkit.modelingtoolkitize(prob)
prob2 = remake(prob, f=ODEFunction(de; jac=true))

sol = solve(prob, TRBDF2())

for Alg in [Rodas5, Rosenbrock23, TRBDF2, KenCarp4]
  @test Array( solve(prob2, Alg(), tstops=sol.t, adaptive=false) ) ≈ Array( solve(prob, Alg(), tstops=sol.t, adaptive=false) ) atol=1e-4
end

## check chunk_size handling in ForwardDiff Jacobians
const chunksize = 1
function rober(du,u,p,t)

  y₁,y₂,y₃ = u
  k₁,k₂,k₃,check = p
  if check && eltype(u) <: ForwardDiff.Dual && ForwardDiff.npartials(u[1])!=chunksize
    @show ForwardDiff.npartials(u[1]), chunksize
    error("chunk_size is not as specifed")
  end


  du[1] = -k₁*y₁+k₃*y₂*y₃
  du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  du[3] =  k₂*y₂^2
  nothing
end
prob1 = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4,true))
sol1 = solve(prob1, TRBDF2(chunk_size=chunksize))
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),(0.04,3e7,1e4,false))
sol = solve(prob, TRBDF2())
@test sol.u[end]==sol1.u[end]
@test length(sol.t)==length(sol1.t)


### DAE Jacobians


function f(out,du,u,p,t)
    out[1] = - 0.04u[1]              + 1e4*u[2]*u[3] - du[1]
    out[2] = + 0.04u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
end

used = Ref(false)
function jac(J, du, u, p, gamma, t)
    used[] = true
    J[1,1] = -0.04 - gamma
    J[1,2] = 1e4*u[3]
    J[1,3] = 1e4*u[2]
    J[2,1] = 0.04
    J[2,2] = -6e7*u[2] - 1e4*u[3] - gamma
    J[2,3] = -1e4*u[2]
    J[3,1] = 1.0
    J[3,2] = 1.0
    J[3,3] = 1.0
end

u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0,100000.0)

differential_vars = [true,true,false]
prob = DAEProblem(DAEFunction(f; jac),du₀,u₀,tspan,differential_vars=differential_vars)
prob2 = DAEProblem(f,du₀,u₀,tspan,differential_vars=differential_vars)

used[] = false
sol = solve(prob, DABDF2())
@test used[]
used[] = false
sol2 = solve(prob2, DABDF2())
@test !used[]
@test iszero(maximum(Array(sol)-Array(sol2)))
