mutable struct AutoSwitch{nAlg,sAlg,tolType}
  count::Int
  nonstiffalg::nAlg
  stiffalg::sAlg
  is_stiffalg::Bool
  maxstiffstep::Int
  maxnonstiffstep::Int
  tol::tolType
end
AutoSwitch(nonstiffalg::nAlg, stiffalg::sAlg; maxstiffstep=15, maxnonstiffstep=15, tol::T=11//10) where {nAlg,sAlg,T} = AutoSwitch(0, nonstiffalg, stiffalg, false, maxstiffstep, maxnonstiffstep, tol)

function is_stiff(eigen_est, dt, alg, tol)
  stiffness = eigen_est*dt/alg_stability_size(alg())
  stiffness < oneunit(stiffness) * tol
end

function (AS::AutoSwitch)(integrator)
  eigen_est, dt = integrator.eigen_est, integrator.dt
  if (AS.count > AS.maxstiffstep && !AS.is_stiffalg)
    AS.is_stiffalg = true
  elseif (AS.count < -AS.maxnonstiffstep && AS.is_stiffalg)
    AS.is_stiffalg = false
  end
  if is_stiff(eigen_est, dt, AS.nonstiffalg, AS.tol)
    AS.count += 1
  else
    AS.count = AS.count > 0 ? -1 : AS.count - 1
  end
  return Int(AS.is_stiffalg) + 1
end

AutoRodas5(alg=Tsit5(); kwargs...) = CompositeAlgorithm((alg, Rodas5()), AutoSwitch(typeof(alg), Rodas5; kwargs...))
