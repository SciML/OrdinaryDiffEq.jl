function initialize!(integrator, cache::BS3ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::BS3ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack a21,a32,a41,a42,a43,c1,c2,btilde1,btilde2,btilde3,btilde4 = cache
  k1 = integrator.fsalfirst
  a1 = dt*a21
  k2 = f(uprev+a1*k1, p, t+c1*dt)
  a2 = dt*a32
  k3 = f(uprev+a2*k2, p, t+c2*dt)
  u = uprev+dt*(a41*k1+a42*k2+a43*k3)
  k4 = f(u, p, t+dt); integrator.fsallast = k4
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator, cache::BS3Cache)
  integrator.kshortsize = 2
  resize!(integrator.k, integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k4
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

@muladd function perform_step!(integrator, cache::BS3Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k2,k3,k4,utilde,tmp,atmp = cache
  @unpack a21,a32,a41,a42,a43,c1,c2,btilde1,btilde2,btilde3,btilde4 = cache.tab
  k1 = cache.fsalfirst
  a1 = dt*a21
  @. tmp = uprev+a1*k1
  f(k2, tmp, p, t+c1*dt)
  a2 = dt*a32
  @. tmp = uprev+a2*k2
  f(k3, tmp, p, t+c2*dt)
  @. u = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, u, p, t+dt)
  if integrator.opts.adaptive
    @. utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator, cache::OwrenZen3ConstantCache)
    integrator.kshortsize = 4
    integrator.k = typeof(integrator.k)(integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:integrator.kshortsize-1
      integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::OwrenZen3ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack a21,a31,a32,a41,a42,a43,c1,c2,btilde1,btilde2,btilde3 = cache
  k1 = integrator.fsalfirst
  a1 = dt*a21
  k2 = f(uprev+a1*k1, p, t+c1*dt)
  tmp =  uprev+ dt*(a31*k1 + a32*k2)
  k3 = f(tmp, p, t+c2*dt)
  u =  uprev+dt*(a41*k1+a42*k2+a43*k3)
  k4 = f(u, p, t+dt); integrator.fsallast = k4
  if integrator.opts.adaptive
    utilde =  dt*(btilde1*k1 + btilde2*k2 + btilde3*k3)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3; integrator.k[4]=k4
  integrator.u = u
end

function initialize!(integrator, cache::OwrenZen3Cache)
    integrator.kshortsize = 4
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
    integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
    integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k4  # setup pointers
    integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

@muladd function perform_step!(integrator, cache::OwrenZen3Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k1,k2,k3,k4,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,c1,c2,btilde1,btilde2,btilde3 = cache.tab
  a1 = dt*a21
  @. tmp = uprev+a1*k1
  f(k2, tmp, p, t+c1*dt)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @. u = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, u, p, t+dt)
  if integrator.opts.adaptive
    @. utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator, cache::OwrenZen4ConstantCache)
  integrator.kshortsize = 6
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::OwrenZen4ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4,btilde1,btilde3,btilde4,btilde5 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(uprev+a*k1, p, t+c1*dt)
  k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
  k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
  k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
  u = uprev+dt*(a61*k1+a63*k3+a64*k4+a65*k5)
  k6 = f(u, p, t+dt); integrator.fsallast = k6
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3; integrator.k[4]=k4;
  integrator.k[5]=k5; integrator.k[6]=k6
  integrator.u = u
end

function initialize!(integrator,cache::OwrenZen4Cache,f=integrator.f)
  integrator.kshortsize = 6
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k6  # setup pointers
  f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
end

#=
@muladd function perform_step!(integrator, cache::OwrenZen4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k1,k2,k3,k4,k5,k6,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4,btilde1,btilde3,btilde4,btilde5 = cache.tab
  a = dt*a21
  @. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, tmp, p, t+c3*dt)
  @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @. u = uprev+dt*(a61*k1+a63*k3+a64*k4+a65*k5)
  f(k6, u, p, t+dt)
  if integrator.opts.adaptive
    @. utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end
=#

@muladd function perform_step!(integrator, cache::OwrenZen4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k1,k2,k3,k4,k5,k6,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4,btilde1,btilde3,btilde4,btilde5 = cache.tab
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(k2, tmp, p, t+c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(k3, tmp, p, t+c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
  f(k4, tmp, p, t+c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(k5, tmp, p, t+c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+dt*(a61*k1[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(k6, u, p, t+dt)
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator, cache::OwrenZen5ConstantCache)
  integrator.kshortsize = 8
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::OwrenZen5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(uprev+a*k1, p, t+c1*dt)
  k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
  k4 = f(uprev+dt*(a41*k1+a42*k2+k3), p, t+c3*dt)
  k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
  k6 = f(uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5), p, t+c5*dt)
  k7 = f(uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6), p, t+c6*dt)
  u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  k8 = f(u, p, t+dt); integrator.fsallast = k8
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3; integrator.k[4]=k4
  integrator.k[5]=k5; integrator.k[6]=k6; integrator.k[7]=k7; integrator.k[8]=k8
  integrator.u = u
end

function initialize!(integrator, cache::OwrenZen5Cache)
  integrator.kshortsize = 8
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;
  integrator.k[7]=cache.k7; integrator.k[8]=cache.k8;
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k8  # setup pointers
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

#=
@muladd function perform_step!(integrator, cache::OwrenZen5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7 = cache.tab
  a = dt*a21
  @. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @. tmp = uprev+dt*(a41*k1+a42*k2+k3)
  f(k4, tmp, p, t+c3*dt)
  @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t+c5*dt)
  @. tmp = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  f(k7, tmp, p, t+c6*dt)
  @. u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  f(k8, u, p, t+dt)
  if integrator.opts.adaptive
    @. utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end
=#

@muladd function perform_step!(integrator, cache::OwrenZen5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7 = cache.tab
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(k2, tmp, p, t+c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(k3, tmp, p, t+c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+k3[i])
  end
  f(k4, tmp, p, t+c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(k5, tmp, p, t+c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(k6, tmp, p, t+c5*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
  end
  f(k7, tmp, p, t+c6*dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
  end
  f(k8, u, p, t+dt)
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator, cache::BS5ConstantCache)
  integrator.kshortsize = 8
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::BS5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(uprev+a*k1, p, t+c1*dt)
  k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
  k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
  k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
  k6 = f(uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5), p, t+c5*dt)
  k7 = f(uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6), p, t+dt)
  u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  integrator.fsallast = f(u, p, t+dt); k8 = integrator.fsallast
  if integrator.opts.adaptive
    uhat   = dt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
    utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
    atmp = calculate_residuals(uhat, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    EEst1 = integrator.opts.internalnorm(atmp)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    EEst2 = integrator.opts.internalnorm(atmp)
    integrator.EEst = max(EEst1,EEst2)
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3;integrator.k[4]=k4;integrator.k[5]=k5;integrator.k[6]=k6;integrator.k[7]=k7;integrator.k[8]=k8
  integrator.u = u
end

function initialize!(integrator, cache::BS5Cache)
  integrator.kshortsize = 8
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;
  integrator.k[7]=cache.k7; integrator.k[8]=cache.k8
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k8  # setup pointers
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

#=
@muladd function perform_step!(integrator, cache::BS5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp = cache
  @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache.tab
  a = dt*a21
  @. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, tmp, p, t+c3*dt)
  @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t+c5*dt)
  @. tmp = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  f(k7, tmp, p, t+dt)
  @. u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  f(k8, u, p, t+dt)
  if integrator.opts.adaptive
    @. utilde = dt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    EEst1 = integrator.opts.internalnorm(atmp)
    @. utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    EEst2 = integrator.opts.internalnorm(atmp)
    integrator.EEst = max(EEst1,EEst2)
  end
end
=#

@muladd function perform_step!(integrator, cache::BS5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp = cache
  @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache.tab
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(k2, tmp, p, t+c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(k3, tmp, p, t+c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
  f(k4, tmp, p, t+c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(k5, tmp, p, t+c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(k6, tmp, p, t+c5*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
  end
  f(k7, tmp, p, t+dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
  end
  f(k8, u, p, t+dt)
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(bhat1*k1[i] + bhat3*k3[i] + bhat4*k4[i] + bhat5*k5[i] + bhat6*k6[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    EEst1 = integrator.opts.internalnorm(atmp)
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i] + btilde8*k8[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    EEst2 = integrator.opts.internalnorm(atmp)
    integrator.EEst = max(EEst1,EEst2)
  end
end

function initialize!(integrator, cache::Tsit5ConstantCache)
  integrator.kshortsize = 7
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::Tsit5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(uprev+a*k1, p, t+c1*dt)
  k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
  k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
  k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
  g6 = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  k6 = f(g6, p, t+dt)
  u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  integrator.fsallast = f(u, p, t+dt); k7 = integrator.fsallast
  if typeof(integrator.alg) <: CompositeAlgorithm
    g7 = u
    # Hairer II, page 22
    integrator.eigen_est = integrator.opts.internalnorm(k7 - k6)/integrator.opts.internalnorm(g7 - g6)
  end
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  integrator.k[1] = k1
  integrator.k[2] = k2
  integrator.k[3] = k3
  integrator.k[4] = k4
  integrator.k[5] = k5
  integrator.k[6] = k6
  integrator.k[7] = k7
  integrator.u = u
end

function initialize!(integrator, cache::Tsit5Cache)
  integrator.kshortsize = 7
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k7 # setup pointers
  resize!(integrator.k, integrator.kshortsize)
  # Setup k pointers
  integrator.k[1] = cache.k1
  integrator.k[2] = cache.k2
  integrator.k[3] = cache.k3
  integrator.k[4] = cache.k4
  integrator.k[5] = cache.k5
  integrator.k[6] = cache.k6
  integrator.k[7] = cache.k7
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

#=
@muladd function perform_step!(integrator, cache::Tsit5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp = cache
  a = dt*a21
  @. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, tmp, p, t+c3*dt)
  @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t+dt)
  @. u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  f(k7, u, p, t+dt)
  if integrator.opts.adaptive
    @. utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end
=#

@muladd function perform_step!(integrator, cache::Tsit5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp = cache
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(k2, tmp, p, t+c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(k3, tmp, p, t+c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
  f(k4, tmp, p, t+c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(k5, tmp, p, t+c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(k6, tmp, p, t+dt)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
  end
  f(k7, u, p, t+dt)
  if typeof(integrator.alg) <: CompositeAlgorithm
    g7 = u
    g6 = tmp
    # Hairer II, page 22
    @. utilde = k7 - k6
    ϱu = integrator.opts.internalnorm(utilde)
    @. utilde = g7 - g6
    ϱd = integrator.opts.internalnorm(utilde)
    integrator.eigen_est = ϱu/ϱd
  end
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde2*k2[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator, cache::DP5ConstantCache)
  integrator.kshortsize = 4
  integrator.k = typeof(integrator.k)(integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  @inbounds for i in eachindex(integrator.k)
    integrator.k[i] = zero(integrator.fsalfirst)
  end
end

@muladd function perform_step!(integrator, cache::DP5ConstantCache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,c1,c2,c3,c4,c5,c6 = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(uprev+a*k1, p, t+c1*dt)
  k3 = f(uprev+dt*(a31*k1+a32*k2), p, t+c2*dt)
  k4 = f(uprev+dt*(a41*k1+a42*k2+a43*k3), p, t+c3*dt)
  k5 = f(uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4), p, t+c4*dt)
  g6 = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  k6 = f(g6, p, t+dt)
  update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
  u = uprev+dt*update
  integrator.fsallast = f(u, p, t+dt); k7 = integrator.fsallast
  if typeof(integrator.alg) <: CompositeAlgorithm
    g7 = u
    # Hairer II, page 22
    integrator.eigen_est = integrator.opts.internalnorm(k7 - k6)/integrator.opts.internalnorm(g7 - g6)
  end
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  integrator.k[1] = update
  bspl = k1 - update
  integrator.k[2] = bspl
  integrator.k[3] = update - k7 - bspl
  integrator.k[4] = d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7
  integrator.u = u
end

function initialize!(integrator, cache::DP5Cache)
  integrator.kshortsize = 4
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = cache.update
  integrator.k[2] = cache.bspl
  integrator.k[3] = cache.dense_tmp3
  integrator.k[4] = cache.dense_tmp4
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k7
  integrator.f(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
end

#=
@muladd function perform_step!(integrator, cache::DP5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,c1,c2,c3,c4,c5,c6 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache.tab
  a = dt*a21
  @. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, tmp, p, t+c3*dt)
  @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t+dt)
  @. update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
  @. u = uprev+dt*update
  f(k7, u, p, t+dt);
  if integrator.opts.adaptive
    @. utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  @. bspl = k1 - update
  @. integrator.k[4] = d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7
  @. integrator.k[3] = update - k7 - bspl
end
=#

@muladd function perform_step!(integrator, cache::DP5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  uidx = eachindex(integrator.uprev)
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,c1,c2,c3,c4,c5,c6 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache.tab
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(k2, tmp, p, t+c1*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(k3, tmp, p, t+c2*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
  f(k4, tmp, p, t+c3*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(k5, tmp, p, t+c4*dt)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(k6, tmp, p, t+dt)
  @tight_loop_macros for i in uidx
    @inbounds update[i] = a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
    @inbounds u[i] = uprev[i]+dt*update[i]
  end
  f(k7, u, p, t+dt)
  if typeof(integrator.alg) <: CompositeAlgorithm
    g6 = tmp
    g7 = u
    # Hairer II, page 22
    @. utilde = k7 - k6
    ϱu = integrator.opts.internalnorm(utilde)
    @. utilde = g7 - g6
    ϱd = integrator.opts.internalnorm(utilde)
    integrator.eigen_est = ϱu/ϱd
  end
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  if integrator.opts.calck
    @tight_loop_macros for i in uidx
      #integrator.k[4] == k5
      @inbounds integrator.k[4][i] = d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k7[i]
      #bspl == k3
      @inbounds bspl[i] = k1[i] - update[i]
      # k6 === integrator.k[3] === k2
      @inbounds integrator.k[3][i] = update[i] - k7[i] - bspl[i]
    end
  end
end
