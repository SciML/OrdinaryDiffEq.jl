mutable struct AutoSwitch{nAlg,sAlg}
  # TODO: Is it better to use UInt8?
  count::Int
  nonstiffalg::nAlg
  stiffalg::sAlg
  is_stiffalg::Bool
  is_successive::Bool
  maxstiffstep::Int
  maxnonstiffstep::Int
end
AutoSwitch(nonstiffalg::nAlg, stiffalg::sAlg; maxstiffstep=15, maxnonstiffstep=15) where {nAlg,sAlg} = AutoSwitch(0, nonstiffalg, stiffalg, false, false, maxstiffstep, maxnonstiffstep)

function is_stiff(eigen_est, dt, alg)
  stiffness = eigen_est*dt/alg_stability_size(alg())
  stiffness < oneunit(stiffness)
end

function (AS::AutoSwitch)(integrator)
  eigen_est, dt = integrator.eigen_est, integrator.dt
  # Switching from a stiff solver to a non-stiff solver is not yet implemented
  if (AS.count > AS.maxstiffstep && AS.is_successive || AS.is_stiffalg)
    AS.count = 0
    AS.is_stiffalg = true
    return 2
  end
  if is_stiff(eigen_est, dt, AS.nonstiffalg)
    AS.count += 1
    if !(AS.is_successive)
      AS.is_successive = true
    end
  else
    AS.is_successive = false
  end
  return 1
end

AutoRodas5(alg=Tsit5(); kwargs...) = CompositeAlgorithm((alg, Rodas5()), AutoSwitch(typeof(alg), Rodas5; kwargs...))
