struct AutoSwitch{Alg}
  count::StiffCount
  alg::Alg
  argnum::UInt8
  is_successive::Bool
end
AutoSwitch(alg::Alg) where Alg = AutoSwitch(zero(StiffCount), alg)

function is_stiff(eigen_est, dt, alg)
  stiffness = eigen_est*dt/alg_stability_size(alg())
  stiffness < oneunit(stiffness)
end

# Need to fix the logic
function (AS::AutoSwitch)(integrator)
  eigen_est, dt = integrator.eigen_est, integrator.dt
  counter = AS.count
  @show integrator.alg
  if (counter.num > 15 || isimplicit(integrator.alg))
    zero!(counter)
    return 2
  end
  is_stiff(eigen_est, dt, AS.alg) && inc!(counter)
  return 1
end

const Tsit5Rodas5 = CompositeAlgorithm((Tsit5(), Rodas5()),
                                       AutoSwitch(Tsit5))
const DP5Rodas5 = CompositeAlgorithm((DP5(), Rodas5()),
                                      AutoSwitch(DP5))
