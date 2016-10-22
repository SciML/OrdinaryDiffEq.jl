using Unitful#, UnitfulPlots

using OrdinaryDiffEq

f = (t,y) -> 0.5*y / 3.0u"s"
u₀ = 1.0u"N"
prob = ODEProblem(f,u₀)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:Midpoint)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:ExplicitRK,adaptive=true)

for alg in OrdinaryDiffEq.DIFFERENTIALEQUATIONSJL_ALGORITHMS
  if !contains(string(alg),"Vectorized") && !contains(string(alg),"Threaded") && alg ∉ OrdinaryDiffEq.DIFFERENTIALEQUATIONSJL_IMPLICITALGS
    sol = solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=alg,adaptive=true)
  end
end

println("Units for Number pass")

TEST_PLOT && plot(sol)

u₀ = [1.0u"N" 2.0u"N"
      3.0u"N" 1.0u"N"]

prob = ODEProblem(f,u₀)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:RK4)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:ExplicitRK)
sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:DP5)

sol =solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=:DP5Threaded)

for alg in OrdinaryDiffEq.DIFFERENTIALEQUATIONSJL_ALGORITHMS
  println(alg)
  if alg ∉ OrdinaryDiffEq.DIFFERENTIALEQUATIONSJL_IMPLICITALGS
    sol = solve(prob::ODEProblem,[0,1],Δt=(1/2^4)u"s",save_timeseries=true,alg=alg,adaptive=true)
  end
end

println("Units for 2D pass")

true
