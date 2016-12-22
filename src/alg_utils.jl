function isfsal(alg::OrdinaryDiffEqAlgorithm)
  if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded || typeof(alg) <: DP8 || typeof(alg) <: BS3 || typeof(alg) <: BS5 || typeof(alg) <: Tsit5 || typeof(alg) <: Vern6 || typeof(alg) <: Rosenbrock23 || typeof(alg) <: Rosenbrock32
    return true
  else
    return false
  end
end

function isspecialdense(alg::OrdinaryDiffEqAlgorithm)
  if typeof(alg) <: DP5 || typeof(alg) <: DP5Threaded || typeof(alg) <: DP8 || typeof(alg) <: BS5 || typeof(alg) <: Tsit5 || typeof(alg) <: Vern6 || typeof(alg) <: Vern7 || typeof(alg) <: Vern8 || typeof(alg) <: Vern9 ||
  typeof(alg) <: Rosenbrock23 || typeof(alg) <: Rosenbrock32
    return true
  else
    return false
  end
end

function isimplicit(alg::OrdinaryDiffEqAlgorithm)
  if typeof(alg) <: ImplicitEuler || typeof(alg) <: Trapezoid
    return true
  else
    return false
  end
end
