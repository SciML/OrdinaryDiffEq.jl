@muladd function perform_step!(integrator, cache::PDIRK44ConstantCache, repeat_step=false)

end

function initialize!(integrator, cache::PDIRK44Cache)
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.fsallast  
end
@muladd function perform_step!(integrator, cache::PDIRK44Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack nlsolver,k1,k2 = cache
  if alg.threading == true
    Threads.@threads for i in 1:2
      indexed_update_W!(integrator, cache, dt, Threads.threadid(), repeat_step)
      _nlsolver = nlsolver[Threads.threadid()]
      _nlsolver.tmp .= uprev
      if Threads.threadid() == 1
        _nlsolver.γ = dt/2
        _nlsolver.c = 1//2
      else
        _nlsolver.γ = 2dt/3
        _nlsolver.c = 2//3
      end
      k1[Threads.threadid()] .= DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    end
    nlsolvefail(nlsolver[1]) && return
    nlsolvefail(nlsolver[2]) && return
    Threads.@threads for i in 1:2
      _nlsolver = nlsolver[Threads.threadid()]
      if Threads.threadid() == 1
        _nlsolver.γ = dt/2
        _nlsolver.c = 1//2
        @.. _nlsolver.tmp .= uprev + dt*(-2.5k1[1]+2.5k1[2])
      else
        _nlsolver.γ = 2dt/3
        _nlsolver.c = 1//3
        @.. _nlsolver.tmp = uprev + dt*((-5//3)k1[1]+(4//3)k1[2])
      end
      k2[Threads.threadid()] .= DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    end
    nlsolvefail(nlsolver[1]) && return
    nlsolvefail(nlsolver[2]) && return
  else
    _nlsolver = nlsolver[1]
    _nlsolver.z .= zero(eltype(u))
    indexed_update_W!(integrator, cache, dt, 1, repeat_step)
    _nlsolver.tmp .= uprev
    _nlsolver.γ = dt/2
    _nlsolver.c = 1//2
    println("here")
    k1[1] = DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    println("here2")
    _nlsolver.tmp .= uprev
    _nlsolver.γ = 2dt/3
    _nlsolver.c = 2//3
    k1[2] = DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    @.. _nlsolver.tmp .= uprev + dt*(-2.5k1[1]+2.5k1[2])
    _nlsolver.γ = dt/2
    _nlsolver.c = 1//2
    k2[1] .= DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    @.. _nlsolver.tmp = uprev + dt*((-5//3)k1[1]+(4//3)k1[2])
    _nlsolver.γ = 2dt/3
    _nlsolver.c = 1//3
    k2[2] .= DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
  end
  @.. u = uprev + dt*(-k1[1] - k2[1] + 1.5k1[2] + 1.5k2[2])
end