function initialize!(integrator, cache::BS3ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

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
  integrator.destats.nf += 3
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::BS3Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k2,k3,k4,utilde,tmp,atmp = cache
  @unpack a21,a32,a41,a42,a43,c1,c2,btilde1,btilde2,btilde3,btilde4 = cache.tab
  # k1 = cache.fsalfirst
  k1 = integrator.fsalfirst
  a1 = dt*a21
  @.. tmp = uprev+a1*k1
  f(k2, tmp, p, t+c1*dt)
  a2 = dt*a32
  @.. tmp = uprev+a2*k2
  f(k3, tmp, p, t+c2*dt)
  @.. u = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, u, p, t+dt)
  integrator.destats.nf += 3
  if integrator.opts.adaptive
    @.. utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end

function initialize!(integrator, cache::OwrenZen3ConstantCache)
    integrator.kshortsize = 4
    integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
    integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.destats.nf += 1

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
  integrator.destats.nf += 3
  if integrator.opts.adaptive
    utilde =  dt*(btilde1*k1 + btilde2*k2 + btilde3*k3)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
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
    integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::OwrenZen3Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k1,k2,k3,k4,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,c1,c2,btilde1,btilde2,btilde3 = cache.tab
  a1 = dt*a21
  @.. tmp = uprev+a1*k1
  f(k2, tmp, p, t+c1*dt)
  @.. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @.. u = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, u, p, t+dt)
  integrator.destats.nf += 3
  if integrator.opts.adaptive
    @.. utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end

function initialize!(integrator, cache::OwrenZen4ConstantCache)
  integrator.kshortsize = 6
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

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
  integrator.destats.nf += 5
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3; integrator.k[4]=k4;
  integrator.k[5]=k5; integrator.k[6]=k6
  integrator.u = u
end

function initialize!(integrator,cache::OwrenZen4Cache)
  integrator.kshortsize = 6
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k6  # setup pointers
  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::OwrenZen4Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k1,k2,k3,k4,k5,k6,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4,btilde1,btilde3,btilde4,btilde5 = cache.tab
  a = dt*a21
  @.. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @.. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @.. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, tmp, p, t+c3*dt)
  @.. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @.. u = uprev+dt*(a61*k1+a63*k3+a64*k4+a65*k5)
  f(k6, u, p, t+dt)
  integrator.destats.nf += 5
  if integrator.opts.adaptive
    @.. utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  return nothing
end

#=
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
  integrator.destats.nf += 5
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end
=#

function initialize!(integrator, cache::OwrenZen5ConstantCache)
  integrator.kshortsize = 8
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

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
  integrator.destats.nf += 7
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::OwrenZen5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7 = cache.tab
  a = dt*a21
  @.. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @.. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @.. tmp = uprev+dt*(a41*k1+a42*k2+k3)
  f(k4, tmp, p, t+c3*dt)
  @.. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @.. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t+c5*dt)
  @.. tmp = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  f(k7, tmp, p, t+c6*dt)
  @.. u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  f(k8, u, p, t+dt)
  integrator.destats.nf += 7
  if integrator.opts.adaptive
    @.. utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  return nothing
end

#=
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
  integrator.destats.nf += 7
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end
=#

function initialize!(integrator, cache::BS5ConstantCache)
  alg = unwrap_alg(integrator, false)
  alg.lazy ? (integrator.kshortsize = 8) : (integrator.kshortsize = 11)
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

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
  integrator.destats.nf += 6
  u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  integrator.fsallast = f(u, p, t+dt); k8 = integrator.fsallast
  integrator.destats.nf += 1
  if integrator.opts.adaptive
    uhat   = dt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
    utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
    atmp = calculate_residuals(uhat, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    EEst1 = integrator.opts.internalnorm(atmp,t)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    EEst2 = integrator.opts.internalnorm(atmp,t)
    integrator.EEst = max(EEst1,EEst2)
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3;integrator.k[4]=k4;integrator.k[5]=k5;integrator.k[6]=k6;integrator.k[7]=k7;integrator.k[8]=k8
  integrator.u = u

  alg = unwrap_alg(integrator, false)
  if !alg.lazy && (integrator.opts.adaptive == false || accept_step_controller(integrator, integrator.opts.controller))
    @unpack c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = cache
    k = integrator.k
    k[9] = f(uprev+dt*(a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8]),p,t+c6*dt)
    k[10] = f(uprev+dt*(a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9]),p,t+c7*dt)
    k[11] = f(uprev+dt*(a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10]),p,t+c8*dt)
    integrator.destats.nf += 3
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::BS5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp = cache
  @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache.tab
  a = dt*a21
  @.. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @.. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @.. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, tmp, p, t+c3*dt)
  @.. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @.. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t+c5*dt)
  @.. tmp = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  f(k7, tmp, p, t+dt)
  @.. u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  f(k8, u, p, t+dt)
  integrator.destats.nf += 7
  if integrator.opts.adaptive
    @.. utilde = dt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    EEst1 = integrator.opts.internalnorm(atmp,t)
    @.. utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    EEst2 = integrator.opts.internalnorm(atmp,t)
    integrator.EEst = max(EEst1,EEst2)
  end
  alg = unwrap_alg(integrator, false)
  if !alg.lazy && (integrator.opts.adaptive == false || accept_step_controller(integrator, integrator.opts.controller))
    k = integrator.k
    @unpack c6,c7,c8,a91,a92,a93,a94,a95,a96,a97,a98,a101,a102,a103,a104,a105,a106,a107,a108,a109,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110 = cache.tab
    @.. tmp = uprev+dt*(a91*k[1]+a92*k[2]+a93*k[3]+a94*k[4]+a95*k[5]+a96*k[6]+a97*k[7]+a98*k[8])
    f(k[9],tmp,p,t+c6*dt)
    @.. tmp = uprev+dt*(a101*k[1]+a102*k[2]+a103*k[3]+a104*k[4]+a105*k[5]+a106*k[6]+a107*k[7]+a108*k[8]+a109*k[9])
    f(k[10],tmp,p,t+c7*dt)
    @.. tmp = uprev+dt*(a111*k[1]+a112*k[2]+a113*k[3]+a114*k[4]+a115*k[5]+a116*k[6]+a117*k[7]+a118*k[8]+a119*k[9]+a1110*k[10])
    f(k[11],tmp,p,t+c8*dt)
    integrator.destats.nf += 3
  end
  return nothing
end

#=
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
  integrator.destats.nf += 7
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(bhat1*k1[i] + bhat3*k3[i] + bhat4*k4[i] + bhat5*k5[i] + bhat6*k6[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    EEst1 = integrator.opts.internalnorm(atmp,t)
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i] + btilde8*k8[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    EEst2 = integrator.opts.internalnorm(atmp,t)
    integrator.EEst = max(EEst1,EEst2)
  end

  alg = unwrap_alg(integrator, false)
  if !alg.lazy && (integrator.opts.adaptive == false || accept_step_controller(integrator, integrator.opts.controller))
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
    integrator.destats.nf += 3
  end
end
=#

function initialize!(integrator, cache::Tsit5ConstantCache)
  integrator.kshortsize = 7
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

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
  integrator.destats.nf += 6
  if typeof(integrator.alg) <: CompositeAlgorithm
    g7 = u
    # Hairer II, page 22
    integrator.eigen_est = integrator.opts.internalnorm(k7 - k6,t)/integrator.opts.internalnorm(g7 - g6,t)
  end
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::Tsit5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp = cache
  a = dt*a21
  @.. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @.. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @.. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, tmp, p, t+c3*dt)
  @.. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @.. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t+dt)
  @.. u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  f(k7, u, p, t+dt)
  integrator.destats.nf += 6
  if integrator.alg isa CompositeAlgorithm
    g7 = u
    g6 = tmp
    # Hairer II, page 22
    @.. utilde = k7 - k6
    ϱu = integrator.opts.internalnorm(utilde,t)
    @.. utilde = g7 - g6
    ϱd = integrator.opts.internalnorm(utilde,t)
    integrator.eigen_est = ϱu/ϱd
  end
  if integrator.opts.adaptive
    @.. utilde = dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  return nothing
end

#=
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
  integrator.destats.nf += 6
  if typeof(integrator.alg) <: CompositeAlgorithm
    g7 = u
    g6 = tmp
    # Hairer II, page 22
    @.. utilde = k7 - k6
    ϱu = integrator.opts.internalnorm(utilde,t)
    @.. utilde = g7 - g6
    ϱd = integrator.opts.internalnorm(utilde,t)
    integrator.eigen_est = ϱu/ϱd
  end
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde2*k2[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
end
=#

function initialize!(integrator, cache::DP5ConstantCache)
  integrator.kshortsize = 4
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

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
  integrator.destats.nf += 6
  if typeof(integrator.alg) <: CompositeAlgorithm
    g7 = u
    # Hairer II, page 22
    integrator.eigen_est = integrator.opts.internalnorm(k7 - k6,t)/integrator.opts.internalnorm(g7 - g6,t)
  end
  if integrator.opts.adaptive
    utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
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
  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator, cache::DP5Cache, repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,btilde1,btilde3,btilde4,btilde5,btilde6,btilde7,c1,c2,c3,c4,c5,c6 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache.tab
  a = dt*a21
  @.. tmp = uprev+a*k1
  f(k2, tmp, p, t+c1*dt)
  @.. tmp = uprev+dt*(a31*k1+a32*k2)
  f(k3, tmp, p, t+c2*dt)
  @.. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(k4, tmp, p, t+c3*dt)
  @.. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(k5, tmp, p, t+c4*dt)
  @.. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(k6, tmp, p, t+dt)
  @.. update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
  @.. u = uprev+dt*update
  f(k7, u, p, t+dt)
  integrator.destats.nf += 6
  if integrator.alg isa CompositeAlgorithm
    g6 = tmp
    g7 = u
    # Hairer II, page 22
    ϱu, ϱd = zero(eltype(k7))^2, zero(eltype(g7))^2
    @inbounds for i in eachindex(k7)
      ϱu += (k7[i] - k6[i])^2
      ϱd += (g7[i] - g6[i])^2
    end
    integrator.eigen_est = sqrt(ϱu/ϱd)*oneunit(t)
  end
  if integrator.opts.adaptive
    @.. utilde = dt*(btilde1*k1 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  if integrator.opts.calck
    #integrator.k[4] == k5
    @.. integrator.k[4] = d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7
    #bspl == k3
    @.. bspl = k1 - update
    # k6 === integrator.k[3] === k2
    @.. integrator.k[3] = update - k7 - bspl
  end
  return nothing
end

#=
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
  integrator.destats.nf += 6
  if typeof(integrator.alg) <: CompositeAlgorithm
    g6 = tmp
    g7 = u
    # Hairer II, page 22
    ϱu, ϱd = zero(eltype(k7))^2, zero(eltype(g7))^2
    @inbounds for i in eachindex(k7)
      ϱu += (k7[i] - k6[i])^2
      ϱd += (g7[i] - g6[i])^2
    end
    integrator.eigen_est = sqrt(ϱu/ϱd)*oneunit(t)
  end
  if integrator.opts.adaptive
    @tight_loop_macros for i in uidx
      @inbounds utilde[i] = dt*(btilde1*k1[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i])
    end
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
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
=#

function initialize!(integrator,cache::KYK2014DGSSPRK_3S2_ConstantCache)
 integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
 integrator.destats.nf += 1
 integrator.kshortsize = 2
 integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)

 # Avoid undefined entries if k is an array of arrays
 integrator.fsallast = zero(integrator.fsalfirst)
 return nothing
end

@muladd function perform_step!(integrator,cache::KYK2014DGSSPRK_3S2_ConstantCache,repeat_step=false)
 @unpack t,dt,uprev,u,f,p = integrator
 @unpack α_10, α_20, α_21, α_30, α_32, β_10, β_21, β_30, β_32, c_1, c_2 = cache
 u_1 =  α_10 * uprev + dt*β_10*integrator.fsalfirst
 u_2 = (
  α_20 * uprev +
  α_21 * u_1 + dt*β_21*f(u_1, p, t + c_1*dt)
 )
 integrator.u = (
  α_30 * uprev + dt*β_30*integrator.fsalfirst +
  α_32 * u_2 + dt*β_32*f(u_2, p, t + c_2*dt)
 )
 integrator.k[1] = integrator.fsalfirst
 integrator.k[2] = f(integrator.u, p, t+dt) # For interpolation, then FSAL'd
 integrator.destats.nf += 3
 integrator.fsallast = integrator.k[2]
 return nothing
end

function initialize!(integrator,cache::KYK2014DGSSPRK_3S2_Cache)
 @unpack k,fsalfirst = cache
 integrator.fsalfirst = fsalfirst
 integrator.fsallast = k
 integrator.kshortsize = 2
 resize!(integrator.k, integrator.kshortsize)
 integrator.k[1] = integrator.fsalfirst
 integrator.k[2] = integrator.fsallast
 integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # FSAL for interpolation
 integrator.destats.nf += 1
 return nothing
end

@muladd function perform_step!(integrator,cache::KYK2014DGSSPRK_3S2_Cache,repeat_step=false)
 @unpack t,dt,uprev,u,f,p = integrator
 @unpack k,fsalfirst,u_1,u_2, kk_1, kk_2 = cache
 @unpack α_10, α_20, α_21, α_30, α_32, β_10, β_21, β_30, β_32, c_1, c_2 = cache.tab

 @.. u_1 =  α_10 * uprev + dt*β_10*integrator.fsalfirst
 f(kk_1, u_1, p, t + c_1*dt)
 @.. u_2 = (
  α_20 * uprev +
  α_21 * u_1 + dt*β_21*kk_1
 )
 f(kk_2, u_2, p, t + c_2*dt)
 @.. integrator.u = (
  α_30 * uprev + dt*β_30*integrator.fsalfirst +
  α_32 * u_2 + dt*β_32*kk_2
 )
 f(integrator.k[2], integrator.u, p, t+dt) # For interpolation, then FSAL'd
 integrator.destats.nf += 3
 return nothing
end


function initialize!(integrator, cache::RKO65ConstantCache)
  integrator.kshortsize = 6
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::RKO65ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p = integrator
  @unpack α21, α31, α41, α51, α32, α42, α52, α62, α43, α53, α63, α54, α64, α65, β2, β3, β4, β5, β6, c1, c2, c3, c4, c5, c6 = cache

  #k1=integrator.fsalfirst #f(uprev,p,t)
  k1 = f(uprev, p, t+c1*dt)
  k2 = f(uprev+α21*dt*k1, p, t+c2*dt)
  k3 = f(uprev+α31*dt*k1+α32*dt*k2, p, t+c3*dt)
  k4 = f(uprev+α41*dt*k1+α42*dt*k2+α43*dt*k3, p, t+c4*dt)
  k5 = f(uprev+α51*dt*k1+α52*dt*k2+α53*dt*k3+α54*dt*k4, p, t+c5*dt)
  k6 = f(uprev+α62*dt*k2+α63*dt*k3+α64*dt*k4+α65*dt*k5, p, t+c6*dt)
  u = uprev+dt*(β2*k2+β3*k3+β4*k4+β5*k5+β6*k6)

  integrator.fsallast = f(u, p, t+dt)  # For interpolation, then FSAL'd

  integrator.destats.nf += 6
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3; integrator.k[4]=k4; integrator.k[5]=k5; integrator.k[6]=k6
  integrator.u = u
end


function initialize!(integrator, cache::RKO65Cache)
  @unpack k,fsalfirst = cache
  integrator.kshortsize = 6
  resize!(integrator.k, integrator.kshortsize)

  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k6  # setup pointers

  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;

  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k6  # setup pointers

  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # Pre-start fsal

  integrator.destats.nf += 1

end



@muladd function perform_step!(integrator, cache::RKO65Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p = integrator
  @unpack tmp, k, k1, k2, k3, k4, k5, k6  = cache
  @unpack α21, α31, α41, α51, α32, α42, α52, α62, α43, α53, α63, α54, α64, α65, β2, β3, β4, β5, β6, c1, c2, c3, c4, c5, c6 = cache.tab
  #println("L221: tmp", tmp)
  f(k1, uprev, p, t+c1*dt)
  @.. tmp = uprev+α21*dt*k1
  #println("L224: tmp/k", tmp, k1)
  f(k2, tmp, p, t+c2*dt)
  @.. tmp = uprev+α31*dt*k1+α32*dt*k2
  f(k3, tmp, p, t+c3*dt)
  @.. tmp = uprev+α41*dt*k1+α42*dt*k2+α43*dt*k3
  f(k4, tmp, p, t+c4*dt)
  @.. tmp = uprev+α51*dt*k1+α52*dt*k2+α53*dt*k3+α54*dt*k4
  f(k5, tmp, p, t+c5*dt)
  @.. tmp = uprev+α62*dt*k2+α63*dt*k3+α64*dt*k4+α65*dt*k5
  f(k6, tmp, p, t+c6*dt)

  @.. u = uprev+dt*(β2*k2+β3*k3+β4*k4+β5*k5+β6*k6)
  #println("L238: tmp/u", tmp, u)
  integrator.destats.nf += 6

  #return nothing
end

function initialize!(integrator, cache::FRK65ConstantCache)
  integrator.kshortsize = 9
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator, cache::FRK65ConstantCache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p = integrator
  @unpack α21, α31, α41, α51, α61, α71, α81, α91, α32, α43, α53, α63, α73, α83, α54, α64, α74, α84, α94, α65, α75, α85, α95, α76, α86, α96, α87, α97, α98, β1, β7, β8, β1tilde, β4tilde, β5tilde, β6tilde, β7tilde, β8tilde, β9tilde, c2, c3, c4, c5, c6, c7, c8, c9, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11 = cache
  alg = unwrap_alg(integrator, true)
  ν = alg.omega*dt
  νsq = ν^2
  β4 = (d1 + νsq*(d2 + νsq*(d3 + νsq*(d4 + νsq*(d5 + νsq*(d6 + +νsq*d7))))))/(1 + νsq*(d8 + νsq*(d9 + νsq*(d10 + νsq*(d11 + νsq*(d12 + +νsq*d13))))))
  β5 = (e1 + νsq*(e2 + νsq*(e3 + νsq*(e4 + νsq*(e5 + νsq*e6)))))/(1 + νsq*(e8 + νsq*(e9 + νsq*(e10 + νsq*e11))))
  β6 = (f1 + νsq*(f2 + νsq*(f3 + νsq*(f4 + νsq*(f5 + νsq*f6)))))/(1 + νsq*(f8 + νsq*(f9 + νsq*(f10 + νsq*f11))))

  k1 = integrator.fsalfirst
  k2 = f(uprev+α21*dt*k1, p, t+c2*dt)
  k3 = f(uprev+α31*dt*k1+α32*dt*k2, p, t+c3*dt)
  k4 = f(uprev+α41*dt*k1+α43*dt*k3, p, t+c4*dt)
  k5 = f(uprev+α51*dt*k1+α53*dt*k3+α54*dt*k4, p, t+c5*dt)
  k6 = f(uprev+α61*dt*k1+α63*dt*k3+α64*dt*k4+α65*dt*k5, p, t+c6*dt)
  k7 = f(uprev+α71*dt*k1+α73*dt*k3+α74*dt*k4+α75*dt*k5+α76*dt*k6, p, t+c7*dt)
  k8 = f(uprev+α81*dt*k1+α83*dt*k3+α84*dt*k4+α85*dt*k5+α86*dt*k6+α87*dt*k7, p, t+c8*dt)
  u = uprev+dt*(β1*k1+β4*k4+β5*k5+β6*k6+β7*k7+β8*k8)
  integrator.fsallast = f(u, p, t+dt)
  k9 = integrator.fsallast
  if integrator.opts.adaptive
    utilde = dt*(β1tilde*k1 + β4tilde*k4 + β5tilde*k5 + β6tilde*k6 + β7tilde*k7 + β8tilde*k8 + β9tilde*k9)
    atmp = calculate_residuals(utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  integrator.destats.nf += 8
  integrator.k[1]=k1
  integrator.k[2]=k2
  integrator.k[3]=k3
  integrator.k[4]=k4
  integrator.k[5]=k5
  integrator.k[6]=k6
  integrator.k[7]=k7
  integrator.k[8]=k8
  integrator.k[9]=k9
  integrator.u = u
end


function initialize!(integrator, cache::FRK65Cache)
  integrator.kshortsize = 9
  integrator.fsalfirst = cache.k1
  integrator.fsallast = cache.k9

  resize!(integrator.k, integrator.kshortsize)

  integrator.k[1]=cache.k1
  integrator.k[2]=cache.k2
  integrator.k[3]=cache.k3
  integrator.k[4]=cache.k4
  integrator.k[5]=cache.k5
  integrator.k[6]=cache.k6
  integrator.k[7]=cache.k7
  integrator.k[8]=cache.k8
  integrator.k[9]=cache.k9

  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # Pre-start fsal
  integrator.destats.nf += 1
end



@muladd function perform_step!(integrator, cache::FRK65Cache, repeat_step=false)
  @unpack t, dt, uprev, u, f, p = integrator
  @unpack tmp, k1, k2, k3, k4, k5, k6, k7, k8, k9, utilde, atmp  = cache
  @unpack α21, α31, α41, α51, α61, α71, α81, α91, α32, α43, α53, α63, α73, α83, α54, α64, α74, α84, α94, α65, α75, α85, α95, α76, α86, α96, α87, α97, α98, β1, β7, β8, β1tilde, β4tilde, β5tilde, β6tilde, β7tilde, β8tilde, β9tilde, c2, c3, c4, c5, c6, c7, c8, c9, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11 = cache.tab
  alg = unwrap_alg(integrator, true)

  ν = alg.omega*dt
  νsq = ν^2
  β4 = (d1 + νsq*(d2 + νsq*(d3 + νsq*(d4 + νsq*(d5 + νsq*(d6 + +νsq*d7))))))/(1 + νsq*(d8 + νsq*(d9 + νsq*(d10 + νsq*(d11 + νsq*(d12 + +νsq*d13))))))
  β5 = (e1 + νsq*(e2 + νsq*(e3 + νsq*(e4 + νsq*(e5 + νsq*e6)))))/(1 + νsq*(e8 + νsq*(e9 + νsq*(e10 + νsq*e11))))
  β6 = (f1 + νsq*(f2 + νsq*(f3 + νsq*(f4 + νsq*(f5 + νsq*f6)))))/(1 + νsq*(f8 + νsq*(f9 + νsq*(f10 + νsq*f11))))

  @.. tmp = uprev+α21*dt*k1
  f(k2, tmp, p, t+c2*dt)
  @.. tmp = uprev+α31*dt*k1+α32*dt*k2
  f(k3, tmp, p, t+c3*dt)
  @.. tmp = uprev+α41*dt*k1+α43*dt*k3
  f(k4, tmp, p, t+c4*dt)
  @.. tmp = uprev+α51*dt*k1+α53*dt*k3+α54*dt*k4
  f(k5, tmp, p, t+c5*dt)
  @.. tmp = uprev+α61*dt*k1+α63*dt*k3+α64*dt*k4+α65*dt*k5
  f(k6, tmp, p, t+c6*dt)
  @.. tmp = uprev+α71*dt*k1+α73*dt*k3+α74*dt*k4+α75*dt*k5+α76*dt*k6
  f(k7, tmp, p, t+c7*dt)
  @.. tmp = uprev+α81*dt*k1+α83*dt*k3+α84*dt*k4+α85*dt*k5+α86*dt*k6+α87*dt*k7
  f(k8, tmp, p, t+c8*dt)
  @.. u = uprev+dt*(β1*k1+β4*k4+β5*k5+β6*k6+β7*k7+β8*k8)
  f(k9, u, p, t+dt)

  integrator.destats.nf += 8
  if integrator.opts.adaptive
    @.. utilde = dt*(β1tilde*k1 + β4tilde*k4 + β5tilde*k5 + β6tilde*k6 + β7tilde*k7 + β8tilde*k8 + β9tilde*k9)
    calculate_residuals!(atmp, utilde, uprev, u, integrator.opts.abstol, integrator.opts.reltol,integrator.opts.internalnorm,t)
    integrator.EEst = integrator.opts.internalnorm(atmp,t)
  end
  return nothing
end

function initialize!(integrator,cache::RKMConstantCache)
  integrator.kshortsize = 6
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
	integrator.fsalfirst = integrator.f(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.destats.nf += 1

	# Avoid undefined entries if k is an array of arrays
	integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::RKMConstantCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack α2, α3, α4, α5, α6, β1, β2, β3, β4, β6, c2, c3, c4, c5, c6 = cache

  #k1 = f(uprev, p, t)
  k1 = integrator.fsalfirst
  k2 = f(uprev+α2*dt*k1, p, t+c2*dt)
  k3 = f(uprev+α3*dt*k2, p, t+c3*dt)
  k4 = f(uprev+α4*dt*k3, p, t+c4*dt)
  k5 = f(uprev+α5*dt*k4, p, t+c5*dt)
  k6 = f(uprev+α6*dt*k5, p, t+c6*dt)
  u = uprev+dt*(β1*k1+β2*k2+β3*k3+β4*k4+β6*k6)

  integrator.fsallast = f(u, p, t+dt) #interpolation then FSAL'd
  integrator.destats.nf += 6
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3; integrator.k[4]=k4; integrator.k[5]=k5; integrator.k[6]=k6
  # integrator.k[1] = integrator.fsalfirst
  # integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::RKMCache)
  @unpack k,fsalfirst = cache
  integrator.kshortsize = 6
  resize!(integrator.k, integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;

  integrator.fsalfirst = cache.k1; integrator.fsallast = zero(integrator.fsalfirst)  # setup pointers

  integrator.f(integrator.fsalfirst,integrator.uprev,integrator.p,integrator.t) # Pre-start fsal

  integrator.destats.nf += 1
end

@muladd function perform_step!(integrator,cache::RKMCache,repeat_step=false)
  @unpack t,dt,uprev,u,f,p = integrator
  @unpack tmp,fsalfirst,k,k1,k2,k3,k4,k5,k6 = cache
  @unpack α2,α3,α4,α5,α6,β1,β2,β3,β4,β6,c2,c3,c4,c5,c6 = cache.tab

  @.. tmp = uprev+α2*dt*k1
  f(k2, tmp, p, t+c2*dt)
  @.. tmp = uprev+α3*dt*k2
  f(k3, tmp, p, t+c3*dt)
  @.. tmp = uprev+α4*dt*k3
  f(k4, tmp, p, t+c4*dt)
  @.. tmp = uprev+α5*dt*k4
  f(k5, tmp, p, t+c5*dt)
  @.. tmp = uprev+α6*dt*k5
  f(k6, tmp, p, t+c6*dt)
  @.. u = uprev+dt*(β1*k1+β2*k2+β3*k3+β4*k4+β6*k6)
  f(integrator.fsallast, u, p, t+dt)
  integrator.destats.nf += 6
  return nothing
end
