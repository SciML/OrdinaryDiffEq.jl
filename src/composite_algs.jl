mutable struct AutoSwitch{nAlg,sAlg,tolType,T}
  count::Int
  nonstiffalg::nAlg
  stiffalg::sAlg
  is_stiffalg::Bool
  maxstiffstep::Int
  maxnonstiffstep::Int
  tol::tolType
  dtfac::T
  stiffalgfirst::Bool
end
AutoSwitch(nonstiffalg::nAlg, stiffalg::sAlg;
           maxstiffstep=15, maxnonstiffstep=15, tol::T=(11//10, 1//10),
           dtfac = 2.0, stiffalgfirst=false) where {nAlg,sAlg,T} =
           AutoSwitch(0, nonstiffalg, stiffalg, stiffalgfirst,
                      maxstiffstep, maxnonstiffstep, tol, dtfac,
                      stiffalgfirst)

function is_stiff(eigen_est, dt, alg, tol)
  stiffness = eigen_est*dt/alg_stability_size(alg)
  stiffness > oneunit(stiffness) * tol[1] || stiffness < oneunit(stiffness) * tol[2]
end

function (AS::AutoSwitch)(integrator)
  eigen_est, dt = integrator.eigen_est, integrator.dt
  if (AS.count > AS.maxstiffstep && !AS.is_stiffalg)
    integrator.dt = dt*AS.dtfac
    AS.is_stiffalg = true
  elseif (AS.count < -AS.maxnonstiffstep && AS.is_stiffalg)
    integrator.dt = dt/AS.dtfac
    AS.is_stiffalg = false
  end
  if is_stiff(eigen_est, dt, AS.nonstiffalg, AS.tol)
    AS.count += 1
  else
    AS.count = AS.count > 0 ? -1 : AS.count - 1
  end
  algnum = AS.stiffalgfirst ? !AS.is_stiffalg : AS.is_stiffalg
  return Int(algnum) + 1
end

function AutoAlgSwitch(nonstiffalg, stiffalg; stiffalgfirst=false, kwargs...)
  AS = AutoSwitch(nonstiffalg, stiffalg; stiffalgfirst=stiffalgfirst, kwargs...)
  stiffalgfirst ? CompositeAlgorithm((stiffalg, nonstiffalg), AS) :
                  CompositeAlgorithm((nonstiffalg, stiffalg), AS)
end

AutoTsit5(alg; kwargs...) = AutoAlgSwitch(Tsit5(), alg; kwargs...)
AutoDP5(alg; kwargs...) = AutoAlgSwitch(DP5(), alg; kwargs...)
