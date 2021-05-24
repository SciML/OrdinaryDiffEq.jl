function initialize!(integrator, cache::PDIRK44ConstantCache) end

@muladd function perform_step!(integrator, cache::PDIRK44ConstantCache, repeat_step=false)
  @unpack dt,uprev,u = integrator
  alg = unwrap_alg(integrator, true)
  @unpack nlsolver, tab = cache
  @unpack γs,cs,α1,α2,b1,b2,b3,b4 = tab

  if isthreaded(alg.threading)
    k2 = Array{typeof(u)}(undef,2)
    k1 = Array{typeof(u)}(undef,2)
    let nlsolver=nlsolver, u=u, uprev=uprev, integrator=integrator, cache=cache, dt=dt, repeat_step=repeat_step,
      k1=k1
      @threaded alg.threading for i in 1:2
        nlsolver[i].z = zero(u)
        nlsolver[i].tmp = uprev
        nlsolver[i].γ = γs[i]
        nlsolver[i].c = cs[i]
        markfirststage!(nlsolver[i])
        k1[i] = nlsolve!(nlsolver[i], integrator, cache, repeat_step)
      end
    end
    nlsolvefail(nlsolver[1]) && return
    nlsolvefail(nlsolver[2]) && return
    let nlsolver=nlsolver, u=u, uprev=uprev, integrator=integrator, cache=cache, dt=dt, repeat_step=repeat_step,
      k1=k1, k2=k2
      @threaded alg.threading for i in 1:2
        nlsolver[i].c = cs[2+i]
        nlsolver[i].z = zero(u)
        nlsolver[i].tmp = uprev + α1[i] * k1[1] + α2[i] * k1[2]
        k2[i] = nlsolve!(nlsolver[i], integrator, cache, repeat_step)
      end
    end
    nlsolvefail(nlsolver[1]) && return
    nlsolvefail(nlsolver[2]) && return
    integrator.u = uprev + b1 * k1[1] + b2 * k2[1] + b3 * k1[2] + b4 * k2[2]
  else
    _nlsolver = nlsolver[1]
    _nlsolver.z = zero(u)
    _nlsolver.tmp = uprev
    _nlsolver.γ = γs[1]
    _nlsolver.c = cs[1]
    markfirststage!(_nlsolver)
    k11 = nlsolve!(_nlsolver, integrator, γs[1]*dt, repeat_step)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z = zero(u)
    _nlsolver.tmp = uprev
    _nlsolver.γ = γs[2]
    _nlsolver.c = cs[2]
    markfirststage!(_nlsolver)
    k12 = nlsolve!(_nlsolver, integrator, γs[2]*dt, repeat_step)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z = zero(u)
    _nlsolver.tmp = uprev + α1[1] * k11 + α2[1] * k12
    _nlsolver.γ = γs[1]
    _nlsolver.c = cs[3]
    markfirststage!(_nlsolver)
    k21 = nlsolve!(_nlsolver, integrator, γs[1]*dt, repeat_step)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z = zero(u)
    _nlsolver.tmp = uprev + α1[2] * k11 + α2[2] * k12
    _nlsolver.γ = γs[2]
    _nlsolver.c = cs[4]
    markfirststage!(_nlsolver)
    k22 = nlsolve!(_nlsolver, integrator, γs[2]*dt, repeat_step)
    nlsolvefail(_nlsolver) && return
    integrator.u = uprev + b1 * k11 + b2 * k21 + b3 * k12 + b4 * k22
  end
end

function initialize!(integrator, cache::PDIRK44Cache) end

@muladd function perform_step!(integrator, cache::PDIRK44Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack nlsolver,k1,k2,tab = cache
  @unpack γs,cs,α1,α2,b1,b2,b3,b4 = tab
  if isthreaded(alg.threading)
    let nlsolver=nlsolver, u=u, uprev=uprev, integrator=integrator, cache=cache, dt=dt, repeat_step=repeat_step,
      k1=k1
      @threaded alg.threading for i in 1:2
        nlsolver[i].z .= zero(eltype(u))
        nlsolver[i].tmp .= uprev
        nlsolver[i].γ = γs[i]
        nlsolver[i].c = cs[i]
        markfirststage!(nlsolver[i])
        k1[i] .= nlsolve!(nlsolver[i], integrator, cache, repeat_step)
      end
    end
    nlsolvefail(nlsolver[1]) && return
    nlsolvefail(nlsolver[2]) && return
    let nlsolver=nlsolver, u=u, uprev=uprev, integrator=integrator, cache=cache, dt=dt, repeat_step=repeat_step,
      k1=k1, k2=k2
      @threaded alg.threading for i in 1:2
        nlsolver[i].c = cs[2+i]
        nlsolver[i].z .= zero(eltype(u))
        @.. nlsolver[i].tmp = uprev + α1[i] * k1[1] + α2[i] * k1[2]
        k2[i] .= nlsolve!(nlsolver[i], integrator, cache, repeat_step)
      end
    end
    nlsolvefail(nlsolver[1]) && return
    nlsolvefail(nlsolver[2]) && return
  else
    _nlsolver = nlsolver[1]
    _nlsolver.z .= zero(eltype(u))
    _nlsolver.tmp .= uprev
    _nlsolver.γ = γs[1]
    _nlsolver.c = cs[1]
    markfirststage!(_nlsolver)
    k1[1] .= nlsolve!(_nlsolver, integrator, cache, repeat_step)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z .= zero(eltype(u))
    _nlsolver.tmp .= uprev
    _nlsolver.γ = γs[2]
    _nlsolver.c = cs[2]
    markfirststage!(_nlsolver)
    k1[2] .= nlsolve!(_nlsolver, integrator, cache, repeat_step)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z .= zero(eltype(u))
    @.. _nlsolver.tmp = uprev + α1[1] * k1[1] + α2[1] * k1[2]
    _nlsolver.γ = γs[1]
    _nlsolver.c = cs[3]
    markfirststage!(_nlsolver)
    k2[1] .= nlsolve!(_nlsolver, integrator, cache, repeat_step)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z .= zero(eltype(u))
    @.. _nlsolver.tmp = uprev + α1[2] * k1[1] + α2[2] * k1[2]
    _nlsolver.γ = γs[2]
    _nlsolver.c = cs[4]
    markfirststage!(_nlsolver)
    k2[2] .= nlsolve!(_nlsolver, integrator, cache, repeat_step)
    nlsolvefail(_nlsolver) && return
  end
  @.. u = uprev + b1 * k1[1] + b2 * k2[1] + b3 * k1[2] + b4 * k2[2]
end
