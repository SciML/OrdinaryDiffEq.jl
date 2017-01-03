using Unitful#, UnitfulPlots
using OrdinaryDiffEq, DiffEqBase
NON_IMPLICIT_ALGS = filter((x)->isleaftype(x) && !OrdinaryDiffEq.isimplicit(x()),union(subtypes(OrdinaryDiffEqAlgorithm),subtypes(OrdinaryDiffEqAdaptiveAlgorithm)))


f = (t,y) -> 0.5*y / 3.0u"s"
u0 = 1.0u"N"
prob = ODEProblem(f,u0,(0.0u"s",1.0u"s"))

sol =solve(prob,Euler(),dt=1u"s"/10)

sol =solve(prob,Midpoint(),dt=1u"s"/10)

sol =solve(prob,ExplicitRK())

for alg in NON_IMPLICIT_ALGS
  if !(alg <: DP5Threaded) && !(alg <: Rosenbrock23) && !(alg <: Rosenbrock32)
    @show alg
    sol = solve(prob,alg(),dt=1u"s"/10)
  end
end

println("Units for Number pass")

u0 = [1.0u"N" 2.0u"N"
      3.0u"N" 1.0u"N"]

prob = ODEProblem(f,u0,(0.0u"s",1.0u"s"))

sol =solve(prob,RK4(),dt=1u"s"/10)

sol =solve(prob,ExplicitRK())
sol =solve(prob,DP5())

sol =solve(prob,DP5Threaded())

for alg in NON_IMPLICIT_ALGS
  if !(alg <: Rosenbrock23) && !(alg <: Rosenbrock32)
    @show alg
    sol = solve(prob,alg(),dt=1u"s"/10)
  end
end

println("Units for 2D pass")
