function initialize!(integrator,cache::KuttaPRK2p5ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.t, integrator) # Pre-start fsal
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::KuttaPRK2p5ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack α21,α31,α32,α41,α42,α43,α5_6 = cache
  @unpack β1,β3,β5,β6,c2,c3,c4,c5_6 = cache

  k1 = f(uprev, t, integrator)
  k2 = f(uprev + dt*α21*k1, t + c2*dt, integrator)
  k3 = f(uprev + dt*(α31*k1 + α32*k2), t + c3*dt, integrator)
  k4 = f(uprev + dt*(α41*k1 + α42*k2 + α43*k3), t + c4*dt, integrator)

  k5_6 = Array{typeof(k1)}(undef, 2)

  if integrator.alg.threading == false
    k5_6[1] = f(uprev + dt*(α5_6[1,1]*k1 + α5_6[1,2]*k2 + α5_6[1,3]*k3 + α5_6[1,4]*k4), t + c5_6[1]*dt, integrator)
    k5_6[2] = f(uprev + dt*(α5_6[2,1]*k1 + α5_6[2,2]*k2 + α5_6[2,3]*k3 + α5_6[2,4]*k4), t + c5_6[2]*dt, integrator)
  else
    let
      Threads.@threads for i in [1,2]
        k5_6[i] = f(uprev + dt*(α5_6[i,1]*k1 + α5_6[i,2]*k2 + α5_6[i,3]*k3 + α5_6[i,4]*k4), t + c5_6[i]*dt, integrator)
      end
    end
  end

  u = uprev + dt*(β1*k1 + β3*k3 + β5*k5_6[1] + β6*k5_6[2])
  k = f(u, t+dt, integrator)

  integrator.fsallast = k # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.k[1] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::KuttaPRK2p5Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.t, integrator) # FSAL for interpolation
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::KuttaPRK2p5Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack k,k1,k2,k3,k4,k5_6,fsalfirst,tmp = cache
  @unpack α21,α31,α32,α41,α42,α43,α5_6 = cache.tab
  @unpack β1,β3,β5,β6,c2,c3,c4,c5_6 = cache.tab

  f( k1, uprev, t, integrator)

  @.. u = uprev + dt*α21*k1
  f( k2, u, t + c2*dt, integrator)

  @.. u = uprev + dt*(α31*k1 + α32*k2)
  f( k3, u, t + c3*dt, integrator)

  @.. u = uprev + dt*(α41*k1 + α42*k2 + α43*k3)
  f( k4, u, t + c4*dt, integrator)

  if integrator.alg.threading == false
    @.. u = uprev + dt*(α5_6[1,1]*k1 + α5_6[1,2]*k2 + α5_6[1,3]*k3 + α5_6[1,4]*k4)
    f( k5_6[1], u, t + c5_6[1]*dt, integrator)

    @.. u = uprev + dt*(α5_6[2,1]*k1 + α5_6[2,2]*k2 + α5_6[2,3]*k3 + α5_6[2,4]*k4)
    f( k5_6[2], u, t + c5_6[2]*dt, integrator)
  else
    tmps = (u, tmp)
    let
      Threads.@threads for i in [1,2]
        @.. tmps[i] = uprev + dt*(α5_6[i,1]*k1 + α5_6[i,2]*k2 + α5_6[i,3]*k3 + α5_6[i,4]*k4)
        f( k5_6[i], tmps[i], t + c5_6[i]*dt, integrator)
      end
    end
  end

  @.. u = uprev + dt*(β1*k1 + β3*k3 + β5*k5_6[1] + β6*k5_6[2])
  f( k, u, t+dt, integrator)
end
