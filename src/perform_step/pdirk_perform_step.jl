@muladd function perform_step!(integrator, cache::PDIRK44ConstantCache, repeat_step=false)

end

function initialize!(integrator, cache::PDIRK44Cache)
  integrator.fsalfirst = similar(cache.k1[1])
  integrator.fsallast = similar(cache.k1[1])  
end
@muladd function perform_step!(integrator, cache::PDIRK44Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p,alg = integrator
  @unpack nlsolver,k1,k2 = cache
  if alg.threading == true
    let nlsolver=nlsolver, u=u, uprev=uprev, integrator=integrator, cache=cache, dt=dt, repeat_step=repeat_step,
      k1=k1
      Threads.@threads for i in 1:2
        nlsolver[i].z .= zero(eltype(u))
        nlsolver[i].tmp .= uprev
        if i == 1
          indexed_update_W!(integrator, cache, dt/2, 1, repeat_step)
          nlsolver[i].γ = 1//2
          nlsolver[i].c = 1//2
        else
          indexed_update_W!(integrator, cache, 2dt/3, 2, repeat_step)
          nlsolver[i].γ = 2//3
          nlsolver[i].c = 2//3
        end
        k1[i] .= DiffEqBase.nlsolve!(nlsolver[i], nlsolver[i].cache, integrator)
      end
    end
    nlsolvefail(nlsolver[1]) && return
    nlsolvefail(nlsolver[2]) && return
    let nlsolver=nlsolver, u=u, uprev=uprev, integrator=integrator, cache=cache, dt=dt, repeat_step=repeat_step,
      k1=k1, k2=k2
      Threads.@threads for i in 1:2
        if i == 1
          nlsolver[i].γ = 1//2
          nlsolver[i].c = 1//2
          nlsolver[i].z .= zero(eltype(u))
          @.. nlsolver[i].tmp = uprev - 2.5 * k1[1] + 2.5 * k1[2]
        else
          nlsolver[i].γ = 2//3
          nlsolver[i].c = 1//3
          nlsolver[i].z .= zero(eltype(u))
          @.. nlsolver[i].tmp = uprev + (-5//3) * k1[1] + (4//3) * k1[2]
        end
        k2[i] .= DiffEqBase.nlsolve!(nlsolver[i], nlsolver[i].cache, integrator)
      end
    end
    nlsolvefail(nlsolver[1]) && return
    nlsolvefail(nlsolver[2]) && return
  else
    _nlsolver = nlsolver[1]
    _nlsolver.z .= zero(eltype(u))
    indexed_update_W!(integrator, cache, dt/2, 1, repeat_step)
    _nlsolver.tmp .= uprev
    _nlsolver.γ = 1//2
    _nlsolver.c = 1//2
    k1[1] .= DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z .= zero(eltype(u))
    indexed_update_W!(integrator, cache, 2dt/3, 1, repeat_step)
    _nlsolver.tmp .= uprev
    _nlsolver.γ = 2//3
    _nlsolver.c = 2//3
    k1[2] .= DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z .= zero(eltype(u))
    indexed_update_W!(integrator, cache, dt/2, 1, repeat_step)
    @.. _nlsolver.tmp .= uprev - 2.5 * k1[1] + 2.5 * k1[2]
    _nlsolver.γ = 1//2
    _nlsolver.c = 1//2
    k2[1] .= DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    nlsolvefail(_nlsolver) && return
    _nlsolver.z .= zero(eltype(u))
    indexed_update_W!(integrator, cache, 2dt/3, 1, repeat_step)
    @.. _nlsolver.tmp = uprev + (-5//3) * k1[1] + (4//3) * k1[2]
    _nlsolver.γ = 2//3
    _nlsolver.c = 1//3
    k2[2] .= DiffEqBase.nlsolve!(_nlsolver, _nlsolver.cache, integrator)
    nlsolvefail(_nlsolver) && return
  end
  @.. u = uprev - k1[1] - k2[1] + 1.5k1[2] + 1.5k2[2]
end