function initialize!(integrator,cache::BS3ConstantCache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::BS3ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4 = cache
  k1 = integrator.fsalfirst
  a1 = dt*a21
  k2 = f(t+c1*dt, @. uprev+a1*k1)
  a2 = dt*a32
  k3 = f(t+c2*dt, @. uprev+a2*k2)
  u = @. uprev+dt*(a41*k1+a42*k2+a43*k3)
  k4 = f(t+dt,u); integrator.fsallast = k4
  if integrator.opts.adaptive
    utilde = @. uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
    tmp = @. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

function initialize!(integrator,cache::BS3Cache,f=integrator.f)
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = cache.fsalfirst  # done by pointers, no copying
  integrator.fsallast = cache.k4
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

@muladd function perform_step!(integrator,cache::BS3Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k2,k3,k4,utilde,tmp,atmp = cache
  @unpack a21,a32,a41,a42,a43,c1,c2,b1,b2,b3,b4 = cache.tab
  k1 = cache.fsalfirst
  a1 = dt*a21
  @. tmp = uprev+a1*k1
  f(t+c1*dt,tmp,k2)
  a2 = dt*a32
  @. tmp = uprev+a2*k2
  f(t+c2*dt,tmp,k3)
  @. u = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(t+dt,u,k4)
  if integrator.opts.adaptive
    @. utilde = uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
    @. atmp = (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator,cache::OwrenZen3ConstantCache,f=integrator.f)
    integrator.kshortsize = 4
    integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
    integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

    # Avoid undefined entries if k is an array of arrays
    integrator.fsallast = zero(integrator.fsalfirst)
    integrator.k[1] = integrator.fsalfirst
    @inbounds for i in 2:integrator.kshortsize-1
      integrator.k[i] = zero(integrator.fsalfirst)
    end
    integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::OwrenZen3ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack a21,a31,a32,a41,a42,a43,c1,c2,b1,b2 = cache
  k1 = integrator.fsalfirst
  a1 = dt*a21
  k2 = f(t+c1*dt, @. uprev+a1*k1)
  tmp = @. uprev+ dt*(a31*k1 + a32*k2)
  k3 = f(t+c2*dt,tmp)
  u = @. uprev+dt*(a41*k1+a42*k2+a43*k3)
  k4 = f(t+dt,u); integrator.fsallast = k4
  if integrator.opts.adaptive
    utilde = @. uprev + dt*(b1*k1 + b2*k2)
    tmp = @. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3; integrator.k[4]=k4
  integrator.u = u
end

function initialize!(integrator,cache::OwrenZen3Cache,f=integrator.f)
    integrator.kshortsize = 4
    integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
    integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
    integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
    integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k4  # setup pointers
    f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

@muladd function perform_step!(integrator,cache::OwrenZen3Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k1,k2,k3,k4,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,c1,c2,b1,b2 = cache.tab
  a1 = dt*a21
  @. tmp = uprev+a1*k1
  f(t+c1*dt,tmp,k2)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(t+c2*dt,tmp,k3)
  @. u = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(t+dt,u,k4)
  if integrator.opts.adaptive
    @. utilde = uprev + dt*(b1*k1 + b2*k2)
    @. atmp = (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator,cache::OwrenZen4ConstantCache,f=integrator.f)
  integrator.kshortsize = 6
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::OwrenZen4ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4,b1,b3,b4 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    uidx = eachindex(uprev)
    tmp = similar(uprev)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+a*k1[i]
    end
    k2 = f(t+c1*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    k3 = f(t+c2*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    k4 = f(t+c3*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    k5 = f(t+c4*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds u[i] = uprev[i]+dt*(a61*k1[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    integrator.fsallast = f(t+dt,u); k6 = integrator.fsallast
    if integrator.opts.adaptive
      atmp = similar(u, typeof(one(recursive_eltype(u))), indices(u))
      @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
        @inbounds utilde   = uprev[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i])
        @inbounds atmp[i] = (utilde-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
      end
      integrator.EEst = integrator.opts.internalnorm(atmp)
    end
  else
    k2 = f(t+c1*dt, uprev+a*k1)
    k3 = f(t+c2*dt, uprev+dt*(a31*k1+a32*k2))
    k4 = f(t+c3*dt, uprev+dt*(a41*k1+a42*k2+a43*k3))
    k5 = f(t+c4*dt, uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
    u = uprev+dt*(a61*k1+a63*k3+a64*k4+a65*k5)
    integrator.fsallast = f(t+dt,u); k6 = integrator.fsallast
    if integrator.opts.adaptive
      utilde = uprev + dt*(b1*k1 + b3*k3 + b4*k4)
      integrator.EEst = integrator.opts.internalnorm(@. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
    end
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3;integrator.k[4]=k4;
  integrator.k[5]=k5;integrator.k[6]=k6
  integrator.u = u
end

function initialize!(integrator,cache::OwrenZen4Cache,f=integrator.f)
  integrator.kshortsize = 6
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k6  # setup pointers
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

@muladd function perform_step!(integrator,cache::OwrenZen4Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k1,k2,k3,k4,k5,k6,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a63,a64,a65,c1,c2,c3,c4,b1,b3,b4 = cache.tab
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(t+c1*dt,tmp,k2)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(t+c2*dt,tmp,k3)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
  f(t+c3*dt,tmp,k4)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(t+c4*dt,tmp,k5)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+dt*(a61*k1[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(t+dt,u,k6)
  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds utilde   = uprev[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i])
      @inbounds atmp[i] = (utilde-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator,cache::OwrenZen5ConstantCache,f=integrator.f)
  integrator.kshortsize = 8
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::OwrenZen5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6,b1,b3,b4,b5,b6 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    uidx = eachindex(uprev)
    tmp = similar(uprev)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+a*k1[i]
    end
    k2 = f(t+c1*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    k3 = f(t+c2*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+k3[i])
    end
    k4 = f(t+c3*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    k5 = f(t+c4*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    k6 = f(t+c5*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
    end
    k7 = f(t+c6*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds u[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
    end
    integrator.fsallast = f(t+dt,u); k8 = integrator.fsallast
    if integrator.opts.adaptive
      atmp = similar(u, typeof(one(recursive_eltype(u))), indices(u))
      @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
        @inbounds utilde   = uprev[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i])
        @inbounds atmp[i] = (utilde-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
      end
      integrator.EEst = integrator.opts.internalnorm(atmp)
    end
  else
    k2 = f(t+c1*dt, uprev+a*k1)
    k3 = f(t+c2*dt, uprev+dt*(a31*k1+a32*k2))
    k4 = f(t+c3*dt, uprev+dt*(a41*k1+a42*k2+k3))
    k5 = f(t+c4*dt, uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
    k6 = f(t+c5*dt, uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
    k7 = f(t+c6*dt, uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6))
    u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
    integrator.fsallast = f(t+dt,u); k8 = integrator.fsallast
    if integrator.opts.adaptive
      utilde = uprev + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6)
      integrator.EEst = integrator.opts.internalnorm(@. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
    end
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3; integrator.k[4]=k4
  integrator.k[5]=k5; integrator.k[6]=k6; integrator.k[7]=k7; integrator.k[8]=k8
  integrator.u = u
end

function initialize!(integrator,cache::OwrenZen5Cache,f=integrator.f)
  integrator.kshortsize = 8
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;
  integrator.k[7]=cache.k7; integrator.k[8]=cache.k8;
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k8  # setup pointers
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

@muladd function perform_step!(integrator,cache::OwrenZen5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,tmp,atmp = cache
  @unpack a21,a31,a32,a41,a42,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,c1,c2,c3,c4,c5,c6,b1,b3,b4,b5,b6 = cache.tab
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(t+c1*dt,tmp,k2)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(t+c2*dt,tmp,k3)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+k3[i])
  end
  f(t+c3*dt,tmp,k4)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(t+c4*dt,tmp,k5)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(t+dt,tmp,k6)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
  end
  f(t+c6*dt,tmp,k7)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
  end
  f(t+dt,u,k8)
  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds utilde   = uprev[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i])
      @inbounds atmp[i] = (utilde-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator,cache::BS5ConstantCache,f=integrator.f)
  integrator.kshortsize = 8
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

#=
@muladd function perform_step!(integrator,cache::BS5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(t+c1*dt, @. uprev+a*k1)
  k3 = f(t+c2*dt, @. uprev+dt*(a31*k1+a32*k2))
  k4 = f(t+c3*dt, @. uprev+dt*(a41*k1+a42*k2+a43*k3))
  k5 = f(t+c4*dt, @. uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
  k6 = f(t+c5*dt, @. uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
  k7 = f(t+dt, @. uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6))
  u = @. uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  integrator.fsallast = f(t+dt,u); k8 = integrator.fsallast
  if integrator.opts.adaptive
    uhat   = @. dt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
    utilde = @. uprev + dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
    EEst1 = integrator.opts.internalnorm(@. uhat/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
    EEst2 = integrator.opts.internalnorm(@. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
    integrator.EEst = max(EEst1,EEst2)
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3;integrator.k[4]=k4;integrator.k[5]=k5;integrator.k[6]=k6;integrator.k[7]=k7;integrator.k[8]=k8
  integrator.u = u
end
=#

@muladd function perform_step!(integrator,cache::BS5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    uidx = eachindex(uprev)
    tmp = similar(uprev)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+a*k1[i]
    end
    k2 = f(t+c1*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    k3 = f(t+c2*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    k4 = f(t+c3*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    k5 = f(t+c4*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    k6 = f(t+c5*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
    end
    k7 = f(t+dt,tmp)
    utmp = similar(u)
    @tight_loop_macros for i in uidx
      @inbounds utmp[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
    end
    u = utmp
    integrator.fsallast = f(t+dt,u); k8 = integrator.fsallast
    if integrator.opts.adaptive
      atmp = similar(u, typeof(one(recursive_eltype(u))), indices(u))
      atmptilde = similar(u, typeof(one(recursive_eltype(u))), indices(u))
      @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
        @inbounds uhat   = dt*(bhat1*k1[i] + bhat3*k3[i] + bhat4*k4[i] + bhat5*k5[i] + bhat6*k6[i])
        @inbounds utilde = uprev[i] + dt*(btilde1*k1[i] + btilde2*k2[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i] + btilde8*k8[i])
        @inbounds atmp[i] = uhat[i]/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
        @inbounds atmptilde[i] = (utilde[i]-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
      end
      EEst1 = integrator.opts.internalnorm(atmp)
      EEst2 = integrator.opts.internalnorm(atmptilde)
      integrator.EEst = max(EEst1,EEst2)
    end
  else
    k2 = f(t+c1*dt, uprev+a*k1)
    k3 = f(t+c2*dt, uprev+dt*(a31*k1+a32*k2))
    k4 = f(t+c3*dt, uprev+dt*(a41*k1+a42*k2+a43*k3))
    k5 = f(t+c4*dt, uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
    k6 = f(t+c5*dt, uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
    k7 = f(t+dt, uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6))
    u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
    integrator.fsallast = f(t+dt,u); k8 = integrator.fsallast
    if integrator.opts.adaptive
      uhat   = dt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
      utilde = uprev + dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
      EEst1 = integrator.opts.internalnorm(@. uhat/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
      EEst2 = integrator.opts.internalnorm(@. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
      integrator.EEst = max(EEst1,EEst2)
    end
  end
  integrator.k[1]=k1; integrator.k[2]=k2; integrator.k[3]=k3;integrator.k[4]=k4;integrator.k[5]=k5;integrator.k[6]=k6;integrator.k[7]=k7;integrator.k[8]=k8
  integrator.u = u
end

function initialize!(integrator,cache::BS5Cache,f=integrator.f)
  integrator.kshortsize = 8
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1]=cache.k1; integrator.k[2]=cache.k2;
  integrator.k[3]=cache.k3; integrator.k[4]=cache.k4;
  integrator.k[5]=cache.k5; integrator.k[6]=cache.k6;
  integrator.k[7]=cache.k7; integrator.k[8]=cache.k8
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k8  # setup pointers
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

#=
@muladd function perform_step!(integrator,cache::BS5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde = cache
  @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache.tab
  a = dt*a21
  @. tmp = uprev+a*k1
  f(t+c1*dt,tmp,k2)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(t+c2*dt,tmp,k3)
  @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(t+c3*dt,tmp,k4)
  @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(t+c4*dt,tmp,k5)
  @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(t+c5*dt,tmp,k6)
  @. tmp = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  f(t+dt,tmp,k7)
  @. u = uprev+dt*(a81*k1+a83*k3+a84*k4+a85*k5+a86*k6+a87*k7)
  f(t+dt,u,k8)
  if integrator.opts.adaptive
    @. uhat   = dt*(bhat1*k1 + bhat3*k3 + bhat4*k4 + bhat5*k5 + bhat6*k6)
    @. utilde = uprev + dt*(btilde1*k1 + btilde2*k2 + btilde3*k3 + btilde4*k4 + btilde5*k5 + btilde6*k6 + btilde7*k7 + btilde8*k8)
    @. atmp = (uhat)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    @. atmptilde = (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    EEst1 = integrator.opts.internalnorm(atmp)
    EEst2 = integrator.opts.internalnorm(atmptilde)
    integrator.EEst = max(EEst1,EEst2)
  end
end
=#

@muladd function perform_step!(integrator,cache::BS5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack k1,k2,k3,k4,k5,k6,k7,k8,utilde,uhat,tmp,atmp,atmptilde = cache
  @unpack c1,c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,a81,a83,a84,a85,a86,a87,bhat1,bhat3,bhat4,bhat5,bhat6,btilde1,btilde2,btilde3,btilde4,btilde5,btilde6,btilde7,btilde8 = cache.tab
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(t+c1*dt,tmp,k2)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(t+c2*dt,tmp,k3)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
  f(t+c3*dt,tmp,k4)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(t+c4*dt,tmp,k5)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(t+c5*dt,tmp,k6)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
  end
  f(t+dt,tmp,k7)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+dt*(a81*k1[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i])
  end
  f(t+dt,u,k8)
  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds uhat[i]   = dt*(bhat1*k1[i] + bhat3*k3[i] + bhat4*k4[i] + bhat5*k5[i] + bhat6*k6[i])
      @inbounds utilde[i] = uprev[i] + dt*(btilde1*k1[i] + btilde2*k2[i] + btilde3*k3[i] + btilde4*k4[i] + btilde5*k5[i] + btilde6*k6[i] + btilde7*k7[i] + btilde8*k8[i])
      @inbounds atmp[i] = uhat[i]/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
      @inbounds atmptilde[i] = (utilde[i]-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    EEst1 = integrator.opts.internalnorm(atmp)
    EEst2 = integrator.opts.internalnorm(atmptilde)
    integrator.EEst = max(EEst1,EEst2)
  end
end

function initialize!(integrator,cache::Tsit5ConstantCache,f=integrator.f)
  integrator.kshortsize = 7
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  @inbounds for i in 2:integrator.kshortsize-1
    integrator.k[i] = zero(integrator.fsalfirst)
  end
  integrator.k[integrator.kshortsize] = integrator.fsallast
end

#=
@muladd function perform_step!(integrator,cache::Tsit5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(t+c1*dt, @. uprev+a*k1)
  k3 = f(t+c2*dt, @. uprev+dt*(a31*k1+a32*k2))
  k4 = f(t+c3*dt, @. uprev+dt*(a41*k1+a42*k2+a43*k3))
  k5 = f(t+c4*dt, @. uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
  k6 = f(t+dt, @. uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
  u = @. uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  integrator.fsallast = f(t+dt,u); k7 = integrator.fsallast
  if integrator.opts.adaptive
    utilde = @. uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
    tmp = @. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
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
=#

@muladd function perform_step!(integrator,cache::Tsit5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    uidx = eachindex(uprev)
    tmp = similar(uprev)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+a*k1[i]
    end
    k2 = f(t+c1*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    k3 = f(t+c2*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    k4 = f(t+c3*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    k5 = f(t+c4*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    k6 = f(t+dt,tmp)
    utmp = similar(u)
    @tight_loop_macros for i in uidx
      @inbounds utmp[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
    end
    u = utmp
    integrator.fsallast = f(t+dt,u); k7 = integrator.fsallast
    if integrator.opts.adaptive
      atmp = similar(u, typeof(one(recursive_eltype(u))), indices(u))
      @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
        @inbounds utilde = uprev[i] + dt*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
        @inbounds atmp[i] = (utilde-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
      end
      integrator.EEst = integrator.opts.internalnorm(atmp)
    end
  else
    k2 = f(t+c1*dt, uprev+a*k1)
    k3 = f(t+c2*dt, uprev+dt*(a31*k1+a32*k2))
    k4 = f(t+c3*dt, uprev+dt*(a41*k1+a42*k2+a43*k3))
    k5 = f(t+c4*dt, uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
    k6 = f(t+dt, uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
    u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    integrator.fsallast = f(t+dt,u); k7 = integrator.fsallast
    if integrator.opts.adaptive
      utilde = uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
      tmp = @. ((utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
      integrator.EEst = integrator.opts.internalnorm(tmp)
    end
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

function initialize!(integrator,cache::Tsit5Cache,f=integrator.f)
  integrator.kshortsize = 7
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k7 # setup pointers
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  # Setup k pointers
  integrator.k[1] = cache.k1
  integrator.k[2] = cache.k2
  integrator.k[3] = cache.k3
  integrator.k[4] = cache.k4
  integrator.k[5] = cache.k5
  integrator.k[6] = cache.k6
  integrator.k[7] = cache.k7
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

#=
@muladd function perform_step!(integrator,cache::Tsit5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp = cache
  a = dt*a21
  @. tmp = uprev+a*k1
  f(t+c1*dt,tmp,k2)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(t+c2*dt,tmp,k3)
  @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(t+c3*dt,tmp,k4)
  @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(t+c4*dt,tmp,k5)
  @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(t+dt,tmp,k6)
  @. u = uprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
  f(t+dt,u,k7)
  if integrator.opts.adaptive
    @. utilde = uprev + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
    @. atmp = (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end
=#

@muladd function perform_step!(integrator,cache::Tsit5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack c1,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76,b1,b2,b3,b4,b5,b6,b7 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,utilde,tmp,atmp = cache
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(t+c1*dt,tmp,k2)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(t+c2*dt,tmp,k3)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
  f(t+c3*dt,tmp,k4)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(t+c4*dt,tmp,k5)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(t+dt,tmp,k6)
  @tight_loop_macros for i in uidx
    @inbounds u[i] = uprev[i]+dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i])
  end
  f(t+dt,u,k7)
  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds utilde[i] = uprev[i] + dt*(b1*k1[i] + b2*k2[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
      @inbounds atmp[i] = (utilde[i]-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

function initialize!(integrator,cache::DP5ConstantCache,f=integrator.f)
  integrator.kshortsize = 4
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.fsalfirst = f(integrator.t,integrator.uprev) # Pre-start fsal

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  @inbounds for i in eachindex(integrator.k)
    integrator.k[i] = zero(integrator.fsalfirst)
  end
end

#=
@muladd function perform_step!(integrator,cache::DP5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  k2 = f(t+c1*dt, @. uprev+a*k1)
  k3 = f(t+c2*dt, @. uprev+dt*(a31*k1+a32*k2))
  k4 = f(t+c3*dt, @. uprev+dt*(a41*k1+a42*k2+a43*k3))
  k5 = f(t+c4*dt, @. uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
  k6 = f(t+dt, @. uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
  update = @. a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
  u = @. uprev+dt*update
  integrator.fsallast = f(t+dt,u); k7 = integrator.fsallast
  if integrator.opts.adaptive
    utilde = @. uprev + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
    tmp = @. (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(tmp)
  end
  integrator.k[1] = update
  bspl = k1 - update
  integrator.k[2] = bspl
  integrator.k[3] = @. update - k7 - bspl
  integrator.k[4] = @. d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7
  integrator.u = u
end
=#

@muladd function perform_step!(integrator,cache::DP5ConstantCache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache
  k1 = integrator.fsalfirst
  a = dt*a21
  if typeof(u) <: AbstractArray && !(typeof(u) <: SArray)
    uidx = eachindex(uprev)
    tmp = similar(uprev)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+a*k1[i]
    end
    k2 = f(t+c1*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
    end
    k3 = f(t+c2*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
    end
    k4 = f(t+c3*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
    end
    k5 = f(t+c4*dt,tmp)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
    end
    k6 = f(t+dt,tmp)
    utmp = similar(u)
    @tight_loop_macros for i in uidx
      @inbounds tmp[i] = a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
      @inbounds utmp[i] = uprev[i]+dt*tmp[i]
    end
    u = utmp
    integrator.fsallast = f(t+dt,u); k7 = integrator.fsallast
    if integrator.opts.adaptive
      atmp = similar(u, typeof(one(recursive_eltype(u))), indices(u))
      @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
        @inbounds utilde = uprev[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
        @inbounds atmp[i] = (utilde-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
      end
      integrator.EEst = integrator.opts.internalnorm(atmp)
    end
    integrator.k[1] = tmp
    k2tmp = similar(integrator.k[2]); k3tmp = similar(integrator.k[3])
    k4tmp = similar(integrator.k[4])
    @tight_loop_macros for i in uidx
      @inbounds k2tmp[i] = k1[i] - tmp[i]
      @inbounds k3tmp[i] = tmp[i] - k7[i] - k2tmp[i]
      @inbounds k4tmp[i] = d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k7[i]
    end
    integrator.k[2] = k2tmp; integrator.k[3] = k3tmp; integrator.k[4] = k4tmp
  else
    k2 = f(t+c1*dt, uprev+a*k1)
    k3 = f(t+c2*dt, uprev+dt*(a31*k1+a32*k2))
    k4 = f(t+c3*dt, uprev+dt*(a41*k1+a42*k2+a43*k3))
    k5 = f(t+c4*dt, uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4))
    k6 = f(t+dt, uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5))
    update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
    u = uprev+dt*update
    integrator.fsallast = f(t+dt,u); k7 = integrator.fsallast
    if integrator.opts.adaptive
      utilde = uprev + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
      tmp = @. ((utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol))
      integrator.EEst = integrator.opts.internalnorm(tmp)
    end
    integrator.k[1] = update
    bspl = k1 - update
    integrator.k[2] = bspl
    integrator.k[3] = update - k7 - bspl
    integrator.k[4] = d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7
  end
  integrator.u = u
end

function initialize!(integrator,cache::DP5Cache,f=integrator.f)
  integrator.kshortsize = 4
  integrator.k = [cache.update,cache.bspl,cache.dense_tmp3,cache.dense_tmp4]
  integrator.fsalfirst = cache.k1; integrator.fsallast = cache.k7
  f(integrator.t,integrator.uprev,integrator.fsalfirst) # Pre-start fsal
end

#=
@muladd function perform_step!(integrator,cache::DP5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache.tab
  a = dt*a21
  @. tmp = uprev+a*k1
  f(t+c1*dt,tmp,k2)
  @. tmp = uprev+dt*(a31*k1+a32*k2)
  f(t+c2*dt,tmp,k3)
  @. tmp = uprev+dt*(a41*k1+a42*k2+a43*k3)
  f(t+c3*dt,tmp,k4)
  @. tmp = uprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
  f(t+c4*dt,tmp,k5)
  @. tmp = uprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
  f(t+dt,tmp,k6)
  @. update = a71*k1+a73*k3+a74*k4+a75*k5+a76*k6
  @. u = uprev+dt*update
  f(t+dt,u,k7);
  if integrator.opts.adaptive
    @. utilde = uprev + dt*(b1*k1 + b3*k3 + b4*k4 + b5*k5 + b6*k6 + b7*k7)
    @. atmp = (utilde-u)/(integrator.opts.abstol+max(abs(uprev),abs(u))*integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  @. bspl = k1 - update
  @. integrator.k[3] = update - k7 - bspl
  @. integrator.k[4] = d1*k1+d3*k3+d4*k4+d5*k5+d6*k6+d7*k7
end
=#

@muladd function perform_step!(integrator,cache::DP5Cache,f=integrator.f)
  @unpack t,dt,uprev,u = integrator
  uidx = eachindex(integrator.uprev)
  @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a73,a74,a75,a76,b1,b3,b4,b5,b6,b7,c1,c2,c3,c4,c5,c6 = cache.tab
  @unpack k1,k2,k3,k4,k5,k6,k7,dense_tmp3,dense_tmp4,update,bspl,utilde,tmp,atmp = cache
  @unpack d1,d3,d4,d5,d6,d7 = cache.tab
  a = dt*a21
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+a*k1[i]
  end
  f(t+c1*dt,tmp,k2)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a31*k1[i]+a32*k2[i])
  end
  f(t+c2*dt,tmp,k3)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a41*k1[i]+a42*k2[i]+a43*k3[i])
  end
  f(t+c3*dt,tmp,k4)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i])
  end
  f(t+c4*dt,tmp,k5)
  @tight_loop_macros for i in uidx
    @inbounds tmp[i] = uprev[i]+dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i])
  end
  f(t+dt,tmp,k6)
  @tight_loop_macros for i in uidx
    @inbounds update[i] = a71*k1[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]
    @inbounds u[i] = uprev[i]+dt*update[i]
  end
  f(t+dt,u,k7)
  if integrator.opts.adaptive
    @tight_loop_macros for (i,atol,rtol) in zip(uidx,Iterators.cycle(integrator.opts.abstol),Iterators.cycle(integrator.opts.reltol))
      @inbounds utilde[i] = uprev[i] + dt*(b1*k1[i] + b3*k3[i] + b4*k4[i] + b5*k5[i] + b6*k6[i] + b7*k7[i])
      @inbounds atmp[i] = (utilde[i]-u[i])/(atol+max(abs(uprev[i]),abs(u[i]))*rtol)
    end
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
  @tight_loop_macros for i in uidx
    @inbounds bspl[i] = k1[i] - update[i]
    @inbounds integrator.k[3][i] = update[i] - k7[i] - bspl[i]
    @inbounds integrator.k[4][i] = d1*k1[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+d7*k7[i]
  end
end
