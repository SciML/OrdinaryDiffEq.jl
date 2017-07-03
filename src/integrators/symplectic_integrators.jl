# http://www.chimica.unipd.it/antonino.polimeno/pubblica/downloads/JChemPhys_101_4062.pdf

@inline function initialize!(integrator,cache::SymplecticEulerCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  # Do the calculation pre
  # So that way FSAL interpolation
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku = integrator.k[1].x[1]
  kdu = integrator.k[2].x[2]
  f[2](integrator.t,uprev,duprev,kdu)
  dt = integrator.dt
  #du = muladd.(integrator.dt,kdu,duprev)
  duidx = eachindex(du)
  @tight_loop_macros for i in duidx
    @inbounds du[i] = muladd(dt,kdu[i],duprev[i])
  end
  f[1](integrator.t,uprev,du,ku)
end

@inline function perform_step!(integrator,cache::SymplecticEulerCache,f=integrator.f)
  @unpack t,dt = integrator
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  duidx = eachindex(du)
  uidx = eachindex(u)
  kuprev = integrator.k[1].x[1]
  ku  = integrator.k[2].x[1]
  kdu = integrator.k[2].x[2]
  #u .= muladd.(dt,kuprev,uprev)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = muladd(dt,kuprev[i],uprev[i])
  end
  # Now actually compute the step
  # Do it at the end for interpolations!
  f[2](integrator.t,uprev,duprev,kdu)
  #du .= muladd.(dt,kdu,duprev)
  @tight_loop_macros for i in duidx
    @inbounds du[i] = muladd(dt,kdu[i],duprev[i])
  end
  f[1](integrator.t,uprev,du,ku)
end

@inline function initialize!(integrator,cache::VelocityVerletCache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::VelocityVerletCache,f=integrator.f)
  @unpack t,dt = integrator
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # x(t+Δt) = x(t) + v(t)*Δt + 1/2*a(t)*Δt^2
  f[2](t,uprev,duprev,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = @muladd uprev[i]+duprev[i]*dt+(1//2*ku[i])*dt^2
  end
  f[2](t+dt,u,duprev,kdu)
  # v(t+Δt) = v(t) + 1/2*(a(t)+a(t+Δt))*Δt
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = @muladd duprev[i] + dt*(1//2*ku[i] + 1//2*kdu[i])
  end
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::Symplectic2Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Symplectic2Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,b1,b2 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::Symplectic3Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Symplectic3Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,b1,b2,b3 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::Symplectic4Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Symplectic4Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,b1,b2,b3,b4 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b4*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a4*kdu[i]
  end
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::Symplectic45Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Symplectic45Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,b1,b2,b3,b4,b5 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b4*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a4*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b5*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  #=
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a5*kdu[i]
  end
  =#
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::Symplectic5Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Symplectic5Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b4*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a4*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b5*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a5*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b6*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a6*kdu[i]
  end
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::Symplectic6Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Symplectic6Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,b1,b2,b3,b4,b5,b6,b7,b8 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b4*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a4*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b5*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a5*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b6*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a6*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b7*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a7*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b8*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  #=
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a8*kdu[i]
  end
  =#
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::Symplectic62Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Symplectic62Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b4*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a4*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b5*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a5*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b6*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a6*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b7*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a7*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b8*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a8*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b9*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a9*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b10*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  #=
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a10*kdu[i]
  end
  =#
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::McAte8Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::McAte8Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b4*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a4*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b5*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a5*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b6*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a6*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b7*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a7*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b8*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a8*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b9*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a9*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b10*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a10*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b11*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a11*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b12*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a12*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b13*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a13*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b14*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a14*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b15*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a15*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b16*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  #=
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a16*kdu[i]
  end
  =#
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::KahanLi8Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::KahanLi8Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b4*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a4*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b5*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a5*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b6*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a6*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b7*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a7*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b8*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a8*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b9*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a9*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b10*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a10*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b11*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a11*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b12*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a12*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b13*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a13*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b14*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a14*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b15*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a15*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b16*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a16*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b17*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a17*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b18*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  #=
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a18*kdu[i]
  end
  =#
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end

@inline function initialize!(integrator,cache::SofSpa10Cache,f=integrator.f)
  integrator.kshortsize = 2
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  uprev,duprev = integrator.uprev.x
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::SofSpa10Cache,f=integrator.f)
  @unpack t,dt = integrator
  @unpack a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,
          a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,
          a35,a36,
          b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,
          b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,
          b35,b36 = cache.tab
  uprev,duprev = integrator.uprev.x
  u,du = integrator.u.x
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  # update position
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] = uprev[i]+dt*b1*duprev[i]
  end
  # update velocity
  f[2](integrator.t,u,duprev,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] = duprev[i] + dt*a1*kdu[i]
  end
  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b2*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a2*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b3*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a3*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b4*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a4*kdu[i]
  end

  # update position & velocity
  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b5*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a5*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b6*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a6*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b7*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a7*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b8*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a8*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b9*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a9*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b10*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a10*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b11*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a11*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b12*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a12*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b13*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a13*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b14*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a14*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b15*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a15*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b16*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a16*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b17*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a17*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b18*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a18*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b19*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a19*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b20*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a20*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b21*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a21*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b22*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a22*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b23*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a23*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b24*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a24*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b25*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a25*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b26*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a26*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b27*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a27*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b28*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a28*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b29*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a29*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b30*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a30*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b31*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a31*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b32*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a32*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b33*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a33*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b34*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a34*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b35*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a35*kdu[i]
  end

  f[1](integrator.t,u,du,ku)
  @tight_loop_macros for i in eachindex(u)
    @inbounds u[i] += dt*b36*ku[i]
  end

  f[2](integrator.t,u,du,kdu)
  #=
  @tight_loop_macros for i in eachindex(du)
    @inbounds du[i] += dt*a30*kdu[i]
  end
  =#
  copy!(integrator.k[1].x[1],integrator.k[2].x[1])
  copy!(integrator.k[1].x[2],integrator.k[2].x[2])
  copy!(integrator.k[2].x[1],du)
  copy!(integrator.k[2].x[2],kdu)
end
