mutable struct AutoSwitch{Alg}
  count::UInt8
  alg::Alg
  algnum::UInt8
  is_successive::Bool
end
AutoSwitch(alg::Alg) where Alg = AutoSwitch(zero(UInt8), alg,
                                            one(UInt8), false)

function is_stiff(eigen_est, dt, alg)
  stiffness = eigen_est*dt/alg_stability_size(alg())
  stiffness < oneunit(stiffness)
end

function (AS::AutoSwitch)(integrator)
  eigen_est, dt = integrator.eigen_est, integrator.dt
  # Here we assume that the first algorithm is always for non-stiff problems
  if (AS.count > 15 && AS.is_successive || AS.algnum == 2)
    AS.count = 0
    AS.algnum= 2
    return 2
  end
  if is_stiff(eigen_est, dt, AS.alg)
    AS.count += 1
    if !(AS.is_successive)
      AS.is_successive = true
    end
  else
    AS.is_successive = false
  end
  return 1
end

AutoRodas5(;alg=Tsit5) = CompositeAlgorithm((alg(), Rodas5()), AutoSwitch(alg))
