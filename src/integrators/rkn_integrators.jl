## y'' = f(t, y, y')
## y(t₀) = y₀; y'(t₀) = y'₀
## kᵢ' = f(t₀+cᵢh, y₀+cᵢhy'₀+h²∑āᵢⱼk'ⱼ, y'₀+h∑aᵢⱼk'ⱼ)
## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
## y'₁ = y'₀ + h∑bᵢk'ᵢ

const NystromDefaultInitialization = Union{Nystrom4Cache,
                                           Nystrom4VelocityIndependentCache,
                                           Nystrom5VelocityIndependentCache,
                                           IRKN3Cache, IRKN4Cache,
                                           DPRKN6Cache}

function initialize!(integrator,cache::NystromDefaultInitialization,f=integrator.f)
  @unpack fsalfirst,k = cache
  uprev,duprev = integrator.uprev.x

  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f.f1(integrator.t,uprev,duprev,integrator.k[2].x[1])
  f.f2(integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@muladd function perform_step!(integrator,cache::Nystrom4Cache,f=integrator.f)
  @unpack t,dt = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  k₁ = fsalfirst
  halfdt = dt/2
  dtsq = dt^2
  eighth_dtsq = dtsq/8
  half_dtsq = dtsq/2
  ttmp = t+halfdt

  ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
  @. ku = uprev + halfdt*duprev + eighth_dtsq*k₁.x[2]
  ## y'₁ = y'₀ + h∑bᵢk'ᵢ
  @. kdu = duprev + halfdt*k₁.x[2]

  f.f2(ttmp,ku,kdu,k₂.x[2])
  @. ku = uprev + halfdt*duprev + eighth_dtsq*k₁.x[2]
  @. kdu = duprev + halfdt*k₂.x[2]

  f.f2(ttmp,ku,kdu,k₃.x[2])
  @. ku = uprev + dt*duprev + half_dtsq*k₃.x[2]
  @. kdu = duprev + dt*k₃.x[2]

  f.f2(t+dt,ku,kdu,k₄.x[2])
  @. u = uprev + (dtsq/6)*(k₁.x[2] + k₂.x[2] + k₃.x[2]) + dt*duprev
  @. du = duprev + (dt/6)*(k₁.x[2] + k₄.x[2] + 2*(k₂.x[2] + k₃.x[2]))

  f.f1(t+dt,u,du,k.x[1])
  f.f2(t+dt,u,du,k.x[2])
end

@muladd function perform_step!(integrator,cache::Nystrom4VelocityIndependentCache,f=integrator.f)
  @unpack t,dt = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  @unpack tmp,fsalfirst,k₂,k₃,k = cache
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  k₁ = fsalfirst
  halfdt = dt/2
  dtsq = dt^2
  eighth_dtsq = dtsq/8
  half_dtsq = dtsq/2
  ttmp = t+halfdt

  ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
  @. ku = uprev + halfdt*duprev + eighth_dtsq*k₁.x[2]

  f.f2(ttmp,ku,du,k₂.x[2])
  @. ku = uprev + dt*duprev + half_dtsq*k₂.x[2]

  f.f2(t+dt,ku,du,k₃.x[2])
  @. u = uprev + (dtsq/6)*(k₁.x[2] + 2*k₂.x[2]) + dt*duprev
  @. du = duprev + (dt/6)*(k₁.x[2] + k₃.x[2] + 4*k₂.x[2])

  f.f1(t+dt,u,du,k.x[1])
  f.f2(t+dt,u,du,k.x[2])
end

@muladd function perform_step!(integrator,cache::IRKN3Cache,f=integrator.f)
  # if there's a discontinuity or the solver is in the first step
  if integrator.iter < 2 && !integrator.u_modified
    perform_step!(integrator,integrator.cache.onestep_cache)
  else
    @unpack t,dt,k,tprev = integrator
    u,du = integrator.u.x
    uprev, duprev  = integrator.uprev.x
    uprev2,duprev2 = integrator.uprev2.x
    uidx = eachindex(integrator.uprev.x[1])
    @unpack tmp,fsalfirst,k₂,k = cache
    @unpack bconst1,bconst2,c1,a21,b1,b2,bbar1,bbar2 = cache.tab
    ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k₁ = fsalfirst

    f.f2(t+c1*dt,    uprev, duprev, k.x[1])
    f.f2(tprev+c1*dt,uprev2,duprev2,k.x[2])
    @tight_loop_macros for i in uidx
      @inbounds ku[i]  = uprev[i]  + dt*(c1*duprev[i]  + dt*a21*k.x[1][i])
      @inbounds kdu[i] = uprev2[i] + dt*(c1*duprev2[i] + dt*a21*k.x[2][i])
    end

    f.f2(t+c1*dt,    ku, duprev, k₂.x[1])
    f.f2(tprev+c1*dt,kdu,duprev2,k₂.x[2])
    @tight_loop_macros for i in uidx
      @inbounds u[i]  = uprev[i] + bconst1*dt*duprev[i] + dt*(bconst2*duprev2[i] + dt*bbar2*(k₂.x[1][i]-k₂.x[2][i]))
      @inbounds du[i] = duprev[i] + dt*(b1*k.x[1][i] + bbar1*k.x[2][i] + b2*(k₂.x[1][i]-k₂.x[2][i]))
    end
    f.f1(t+dt,u,du,k.x[1])
    f.f2(t+dt,u,du,k.x[2])
  end # end if
end

@muladd function perform_step!(integrator,cache::IRKN4Cache,f=integrator.f)
  # if there's a discontinuity or the solver is in the first step
  if integrator.iter < 2 && !integrator.u_modified
    perform_step!(integrator,integrator.cache.onestep_cache)
  else
    @unpack t,dt,k,tprev = integrator
    u,du = integrator.u.x
    uprev, duprev  = integrator.uprev.x
    uprev2,duprev2 = integrator.uprev2.x
    uidx = eachindex(integrator.uprev.x[1])
    @unpack tmp,fsalfirst,k₂,k₃,k = cache
    @unpack bconst1,bconst2,c1,c2,a21,a32,b1,b2,b3,bbar1,bbar2,bbar3 = cache.tab
    ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
    k₁ = fsalfirst
    dtsq = dt^2

    f.f2(t+c1*dt,    uprev, duprev, k.x[1])
    f.f2(tprev+c1*dt,uprev2,duprev2,k.x[2])
    @tight_loop_macros for i in uidx
      @inbounds ku[i]  = uprev[i]  + dt*(c1*duprev[i]  + dt*a21*k.x[1][i])
      @inbounds kdu[i]  = uprev2[i]  + dt*(c1*duprev2[i]  + dt*a21*k.x[2][i])
    end

    f.f2(t+c1*dt,    ku, duprev, k₂.x[1])
    f.f2(tprev+c1*dt,kdu,duprev2,k₂.x[2])
    @tight_loop_macros for i in uidx
      @inbounds ku[i]  = uprev[i]  + dt*(c2*duprev[i]  + dt*a32*k₂.x[1][i])
      @inbounds kdu[i]  = uprev2[i]  + dt*(c2*duprev2[i]  + dt*a32*k₂.x[2][i])
    end

    f.f2(t+c2*dt,    ku, duprev, k₃.x[1])
    f.f2(tprev+c2*dt,kdu,duprev2,k₃.x[2])
    @tight_loop_macros for i in uidx
      @inbounds u[i]  = uprev[i] + dt*bconst1*duprev[i] + dt*(bconst2*duprev2[i] + dt*(bbar2*(k₂.x[1][i]-k₂.x[2][i]) + bbar3*(k₃.x[1][i]-k₃.x[2][i])))
      @inbounds du[i] = duprev[i] + dt*(b1*k.x[1][i] + bbar1*k.x[2][i] + b2*(k₂.x[1][i]-k₂.x[2][i]) + b3*(k₃.x[1][i]-k₃.x[2][i]))
    end
    f.f1(t+dt,u,du,k.x[1])
    f.f2(t+dt,u,du,k.x[2])
  end # end if
end

#=
@muladd function perform_step!(integrator,cache::Nystrom5VelocityIndependentCache,f=integrator.f)
  @unpack t,dt = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  k₁ = fsalfirst
  dtsq = dt^2

  f.f2(t+1//5*dt,uprev,duprev,k₁.x[2])
  @. ku = uprev + (1//5*dt)*duprev + (1//50*dtsq)*k₁.x[2]

  f.f2(t+1//5*dt,ku,du,k₂.x[2])
  @. ku = uprev + (2//3*dt)*duprev + (-1//27*dtsq)*k₁.x[2] + (7//27*dtsq)*k₂.x[2]

  f.f2(t+2//3*dt,ku,du,k₃.x[2])
  @. ku = uprev + dt*duprev + (3//10*dtsq)*k₁.x[2] + (-2//35*dtsq)*k₂.x[2] + (9//35*dtsq)*k₃.x[2]

  f.f2(t+dt,ku,du,k₄.x[2])
  @. u  = uprev + dt*duprev + (14//336*dtsq)*k₁.x[2] + (100//336*dtsq)*k₂.x[2] + (54//336*dtsq)*k₃.x[2]
  @. du = duprev[i] + (14//336*dt)*k₁.x[2][i] + (125//336*dt)*k₂.x[2][i] + (162//336*dt)*k₃.x[2][i] + (35//336*dt)*k₄.x[2][i]

  f.f1(t+dt,u,du,k.x[1])
  f.f2(t+dt,u,du,k.x[2])
end
=#

@muladd function perform_step!(integrator,cache::Nystrom5VelocityIndependentCache,f=integrator.f)
  @unpack t,dt = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  uidx = eachindex(integrator.uprev.x[1])
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  @unpack c1, c2, a21, a31, a32, a41, a42, a43, bbar1, bbar2, bbar3, b1, b2, b3, b4 = cache.tab
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  k₁ = fsalfirst

  @. ku = uprev + dt*(c1*duprev + dt*a21*k₁.x[2])

  f.f2(t+c1*dt,ku,du,k₂.x[2])
  @. ku = uprev + dt*(c2*duprev + dt*(a31*k₁.x[2] + a32*k₂.x[2]))

  f.f2(t+c2*dt,ku,du,k₃.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(duprev[i] + dt*(a41*k₁.x[2][i] + a42*k₂.x[2][i] + a43*k₃.x[2][i]))
  end

  f.f2(t+dt,ku,du,k₄.x[2])
  @tight_loop_macros for i in uidx
    @inbounds u[i]  = uprev[i] + dt*(duprev[i] + dt*(bbar1*k₁.x[2][i] + bbar2*k₂.x[2][i] + bbar3*k₃.x[2][i]))
    @inbounds du[i] = duprev[i] + dt*(b1*k₁.x[2][i] + b2*k₂.x[2][i] + b3*k₃.x[2][i] + b4*k₄.x[2][i])
  end
  f.f1(t+dt,u,du,k.x[1])
  f.f2(t+dt,u,du,k.x[2])
end

@muladd function perform_step!(integrator,cache::DPRKN6Cache,f=integrator.f)
  @unpack t,dt = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  @unpack tmp,atmp,fsalfirst,k2,k3,k4,k5,k6,k,utilde = cache
  @unpack c1, c2, c3, c4, c5, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a63, a64, a65, b1, b3, b4, b5, bp1, bp3, bp4, bp5, bp6, btilde1, btilde2, btilde3, btilde4, btilde5, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6 = cache.tab
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  uidx = eachindex(integrator.uprev.x[2])
  k1 = fsalfirst

  @. ku = uprev + dt*(c1*duprev + dt*a21*k1.x[2])

  f.f2(t+dt*c1,ku,du,k2.x[2])
  @. ku = uprev + dt*(c2*duprev + dt*(a31*k1.x[2] + a32*k2.x[2]))

  f.f2(t+dt*c2,ku,du,k3.x[2])
  @. ku = uprev + dt*(c3*duprev + dt*(a41*k1.x[2] + a42*k2.x[2] + a43*k3.x[2]))

  f.f2(t+dt*c3,ku,du,k4.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c4*duprev[i] + dt*(a51*k1.x[2][i] + a52*k2.x[2][i] + a53*k3.x[2][i] + a54*k4.x[2][i]))
  end

  f.f2(t+dt*c4,ku,du,k5.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c5*duprev[i] + dt*(a61*k1.x[2][i]                  + a63*k3.x[2][i] + a64*k4.x[2][i] + a65*k5.x[2][i])) # no a62
  end

  f.f2(t+dt*c5,ku,du,k6.x[2])
  @tight_loop_macros for i in uidx
    @inbounds u[i]  = uprev[i] + dt*(duprev[i] + dt*(b1 *k1.x[2][i] + b3 *k3.x[2][i] + b4 *k4.x[2][i] + b5 *k5.x[2][i])) # b1 -- b5, no b2
    @inbounds du[i] = duprev[i]                + dt*(bp1*k1.x[2][i] + bp3*k3.x[2][i] + bp4*k4.x[2][i] + bp5*k5.x[2][i] + bp6*k6.x[2][i]) # bp1 -- bp6, no bp2
  end
  #=
  @tight_loop_macros for i in uidx
    @inbounds u[i]  = uprev[i] + dt*(duprev[i] + dt*(bhat1*k1.x[2][i] + bhat2*k2.x[2][i] + bhat3*k3.x[2][i]))
    @inbounds du[i] = duprev[i]+ dt*(bphat1*k1.x[2][i] + bphat3*k3.x[2][i] + bphat4*k4.x[2][i] + bphat5*k5.x[2][i] + bphat6*k6.x[2][i])
  end
  =#

  f.f1(t+dt,u,du,k.x[1])
  f.f2(t+dt,u,du,k.x[2])
  if integrator.opts.adaptive
    uhat, duhat = utilde.x
    dtsq = dt^2
    @tight_loop_macros for i in uidx
      @inbounds uhat[i]  = dtsq*(btilde1*k1.x[2][i] + btilde2*k2.x[2][i] + btilde3*k3.x[2][i] + btilde4*k4.x[2][i] + btilde5*k5.x[2][i])
      @inbounds duhat[i] = dt*(bptilde1*k1.x[2][i] + bptilde3*k3.x[2][i] + bptilde4*k4.x[2][i] + bptilde5*k5.x[2][i] + bptilde6*k6.x[2][i])
    end
    calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

