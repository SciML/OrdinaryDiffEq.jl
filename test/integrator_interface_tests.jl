using DiffEqBase: set_t!, set_u!, set_ut!
using OrdinaryDiffEq

# set_X!(integrator, integrator.X) should not change the result.
@testset "Trivial $setter ($alg)" for alg in [RK4, Trapezoid],
                                   setter in [set_t!, set_u!, set_ut!]
  f = (u,p,t) -> (2u)
  prob = ODEProblem(f,1/2,(0.0,2.0))
  t_half = 1.0

  integrator1 = init(prob, alg())
  integrator2 = init(prob, alg())

  step!(integrator1, t_half)  # "inexact" stepping w/o tstops
  if setter === set_t!
    set_t!(integrator1, integrator1.t)
  elseif setter === set_u!
    set_u!(integrator1, integrator1.u)
  else
    set_ut!(integrator1, integrator1.u, integrator1.t)
  end
  solve!(integrator1)

  solve!(integrator2)

  @test integrator1.t == integrator2.t
  if setter in [set_u!, set_ut!] || alg === Trapezoid
    rtol = integrator1.opts.reltol
    atol = integrator1.opts.abstol
    if alg === Trapezoid
      rtol *= 100
    end
    @test integrator1.u ≈ integrator2.u  rtol=rtol atol=atol
  else
    @test integrator1.u == integrator2.u
  end
end

@testset "Resolve with $setter ($alg)" for alg in [RK4, Trapezoid],
                                        setter in [set_t!, set_ut!]
  f = (u,p,t) -> (2u * cos(2π * t))
  integrator1 = init(ODEProblem(f,1/2,(0.0,1.0)), alg())
  integrator2 = init(ODEProblem(f,1/2,(0.0,2.0)), alg())

  solve!(integrator1)
  if setter === set_t!
    set_t!(integrator1, 0)
  else
    set_ut!(integrator1, integrator1.u, 0)
  end
  solve!(integrator1)

  solve!(integrator2)

  rtol = integrator1.opts.reltol
  atol = integrator1.opts.abstol
  if alg === Trapezoid
    rtol *= 100
  end
  @test integrator1.u ≈ integrator2.u  rtol=rtol atol=atol
end
