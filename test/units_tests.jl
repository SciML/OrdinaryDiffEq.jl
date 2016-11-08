using Unitful#, UnitfulPlots

NON_IMPLICIT_ALGS = filter((x)->isleaftype(x) && !OrdinaryDiffEq.isimplicit(x()),union(subtypes(OrdinaryDiffEqAlgorithm),subtypes(OrdinaryDiffEqAdaptiveAlgorithm)))
using OrdinaryDiffEq

f = (t,y) -> 0.5*y / 3.0u"s"
u0 = 1.0u"N"
prob = ODEProblem(f,u0,[0;1])

sol =solve(prob,Midpoint(),dt=(1/2^4)u"s",save_timeseries=true)

sol =solve(prob,ExplicitRK(),dt=(1/2^4)u"s",save_timeseries=true,adaptive=true)

for alg in NON_IMPLICIT_ALGS
  if !(alg <: DP5Threaded) && !(alg <: Rosenbrock23) && !(alg <: Rosenbrock32)
    sol = solve(prob,alg(),dt=(1/2^4)u"s",save_timeseries=true,adaptive=true)
  end
end

println("Units for Number pass")

TEST_PLOT && plot(sol)

u0 = [1.0u"N" 2.0u"N"
      3.0u"N" 1.0u"N"]

prob = ODEProblem(f,u0,[0;1])

sol =solve(prob,RK4(),dt=(1/2^4)u"s",save_timeseries=true)

sol =solve(prob,ExplicitRK(),dt=(1/2^4)u"s",save_timeseries=true)
sol =solve(prob,DP5(),dt=(1/2^4)u"s",save_timeseries=true)

sol =solve(prob,DP5Threaded(),dt=(1/2^4)u"s",save_timeseries=true)

for alg in NON_IMPLICIT_ALGS
  if !(alg <: DP5Threaded) && !(alg <: Rosenbrock23) && !(alg <: Rosenbrock32)
    sol = solve(prob,alg(),dt=(1/2^4)u"s",save_timeseries=true,adaptive=true)
  end
end

println("Units for 2D pass")

true
