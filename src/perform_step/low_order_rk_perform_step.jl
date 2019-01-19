function initialize!(integrator, cache::BS3ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
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
  # k1 = cache.fsalfirst
  k1 = integrator.fsalfirst
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
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
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
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
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
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
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
  alg = unwrap_alg(integrator, false)
  alg.lazy ? (integrator.kshortsize = 8) : (integrator.kshortsize = 11)
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:7
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast


  if !alg.lazy
    @inbounds for i in 9:11
      integrator.k[i] = zero(integrator.fsalfirst)
    end
  end
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

  alg = unwrap_alg(integrator, false)
  if !alg.lazy && (integrator.opts.adaptive == false || integrator.EEst <= 1.0)
    @unpack c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = cache
    k = integrator.k
    k[9] = f(uprev+dt*(a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8]),p,t+c6*dt)
    k[10] = f(uprev+dt*(a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9]),p,t+c7*dt)
    k[11] = f(uprev+dt*(a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10]),p,t+c8*dt)
  end
end

function initialize!(integrator, cache::BS5Cache)
  alg = unwrap_alg(integrator, false)
  alg.lazy ? (integrator.kshortsize = 8) : (integrator.kshortsize = 11)
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;
  integrator.k[7]=cache.k7; integrator.k[8]=cache.k8

  if !alg.lazy
    integrator.k[9]= similar(cache.k1)
    integrator.k[10]= similar(cache.k1)
    integrator.k[11]= similar(cache.k1)
  end

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

  alg = unwrap_alg(integrator, false)
  if !alg.lazy && (integrator.opts.adaptive == false || integrator.EEst <= 1.0)
    k = integrator.k
    @unpack c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = cache.tab
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a91*k[1][i]+a92*k[2][i]+a93*k[3][i]+a94*k[4][i]+a95*k[5][i]+a96*k[6][i]+a97*k[7][i]+a98*k[8][i])
    end
    f(k[9],tmp,p,t+c6*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a101*k[1][i]+a102*k[2][i]+a103*k[3][i]+a104*k[4][i]+a105*k[5][i]+a106*k[6][i]+a107*k[7][i]+a108*k[8][i]+a109*k[9][i])
    end
    f(k[10],tmp,p,t+c7*dt)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a111*k[1][i]+a112*k[2][i]+a113*k[3][i]+a114*k[4][i]+a115*k[5][i]+a116*k[6][i]+a117*k[7][i]+a118*k[8][i]+a119*k[9][i]+a1110*k[10][i])
    end
    f(k[11],tmp,p,t+c8*dt)
  end
end

function initialize!(integrator, cache::Tsit5ConstantCache)
  integrator.kshortsize = 7
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
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
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
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
    ϱu, ϱd = zero(eltype(k7)), zero(eltype(k7))
    @inbounds for i in eachindex(k7)
      ϱu += (k7[i] - k6[i])^2
      ϱd += (g7[i] - g6[i])^2
    end
    integrator.eigen_est = sqrt(ϱu/ϱd)
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

function initialize!(integrator,cache::LDDRK46ConstantCache)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.kshortsize = 1
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
end

@muladd function perform_step!(integrator,cache::LDDRK46ConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α1,α2,α3,α4,α5,β1,β2,β3,β4,β5,β6,c2,c3,c4,c5,c6 = cache

  # u0
  u   = uprev
  tmp = dt*integrator.fsalfirst
  #u1
  u   = u + β1*tmp
  tmp = dt*f(u+α1*tmp, p, t+c2*dt)
  #u2
  u   = u + β2*tmp
  tmp = dt*f(u+α2*tmp, p, t+c3*dt)
  #u3
  u   = u + β3*tmp
  tmp = dt*f(u+α3*tmp, p, t+c4*dt)
  #u4
  u   = u + β4*tmp
  tmp = dt*f(u+α4*tmp, p, t+c5*dt)
  #u5
  u   = u + β5*tmp
  tmp = dt*f(u+α5*tmp, p, t+c6*dt)
  #u6
  u   = u + β6*tmp

  integrator.fsallast = f(u, p, t+dt) # For interpolation, then FSAL'd
  integrator.k[1] = integrator.fsalfirst
  integrator.u = u
end

function initialize!(integrator,cache::LDDRK46Cache)
  @unpack k,fsalfirst = cache
  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 1
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
end

@muladd function perform_step!(integrator,cache::LDDRK46Cache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k,fsalfirst,tmp = cache
  @unpack α1,α2,α3,α4,α5,β1,β2,β3,β4,β5,β6,c2,c3,c4,c5,c6 = cache.tab

  #u0
  @. u   = uprev
  @. tmp = dt*fsalfirst
  #u1
  @. u   = u + β1*tmp
  f(k, u+α1*tmp, p, t+c2*dt)
  @. tmp = dt*k
  #u2
  @. u   = u + β2*tmp
  f(k, u+α2*tmp, p, t+c3*dt)
  @. tmp = dt*k
  #u3
  @. u   = u + β3*tmp
  f(k, u+α3*tmp, p, t+c4*dt)
  @. tmp = dt*k
  #u4
  @. u   = u + β4*tmp
  f(k, u+α4*tmp, p, t+c5*dt)
  @. tmp = dt*k
  #u5
  @. u   = u + β5*tmp
  f(k, u+α5*tmp, p, t+c6*dt)
  @. tmp = dt*k
  #u6
  @. u   = u + β6*tmp

  f( k,  u, p, t+dt)
end
