## y'' = f(t, y, y')
## y(t₀) = y₀; y'(t₀) = y'₀
## kᵢ' = f(t₀+cᵢh, y₀+cᵢhy'₀+h²∑āᵢⱼk'ⱼ, y'₀+h∑aᵢⱼk'ⱼ)
## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
## y'₁ = y'₀ + h∑bᵢk'ᵢ

const NystromDefaultInitialization = Union{Nystrom4Cache,
                                           Nystrom4VelocityIndependentCache,
                                           Nystrom5VelocityIndependentCache,
                                           IRKN3Cache, IRKN4Cache,
                                           DPRKN6Cache, DPRKN8Cache,
                                           DPRKN12Cache}

function initialize!(integrator,cache::NystromDefaultInitialization)
  @unpack fsalfirst,k = cache
  uprev,duprev = integrator.uprev.x

  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.f.f1(integrator.t,uprev,duprev,integrator.k[2].x[1])
  integrator.f.f2(integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@muladd function perform_step!(integrator,cache::Nystrom4Cache,repeat_step=false)
  @unpack t,dt,f = integrator
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

@muladd function perform_step!(integrator,cache::Nystrom4VelocityIndependentCache,repeat_step=false)
  @unpack t,dt,f = integrator
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

@muladd function perform_step!(integrator,cache::IRKN3Cache,repeat_step=false)
  # if there's a discontinuity or the solver is in the first step
  if integrator.iter < 2 && !integrator.u_modified
    perform_step!(integrator,integrator.cache.onestep_cache)
  else
    @unpack t,dt,k,tprev,f = integrator
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

@muladd function perform_step!(integrator,cache::IRKN4Cache,repeat_step=false)
  # if there's a discontinuity or the solver is in the first step
  if integrator.iter < 2 && !integrator.u_modified
    perform_step!(integrator,integrator.cache.onestep_cache)
  else
    @unpack t,dt,k,tprev,f = integrator
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
@muladd function perform_step!(integrator,cache::Nystrom5VelocityIndependentCache,repeat_step=false)
  @unpack t,dt,f = integrator
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

@muladd function perform_step!(integrator,cache::Nystrom5VelocityIndependentCache,repeat_step=false)
  @unpack t,dt,f = integrator
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

function initialize!(integrator,cache::DPRKN6Cache,repeat_step=false)
  @unpack fsalfirst,k = cache
  uprev,duprev = integrator.uprev.x

  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 6
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = cache.fsalfirst
  integrator.k[2] = cache.k2
  integrator.k[3] = cache.k3
  integrator.k[4] = cache.k4
  integrator.k[5] = cache.k5
  integrator.k[6] = cache.k6
  integrator.f.f1(integrator.t,uprev,duprev,integrator.fsallast.x[1])
  integrator.f.f2(integrator.t,uprev,duprev,integrator.fsallast.x[2])
end

@muladd function perform_step!(integrator,cache::DPRKN6Cache,repeat_step=false)
  @unpack t,dt,f = integrator
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

@muladd function perform_step!(integrator,cache::DPRKN8Cache,repeat_step=false)
  @unpack t,dt,f = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  @unpack tmp,atmp,fsalfirst,k2,k3,k4,k5,k6,k7,k8,k9,k,utilde = cache
  @unpack c1, c2, c3, c4, c5, c6, c7, c8, a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65, a71, a72, a73, a74, a75, a76, a81, a82, a83, a84, a85, a86, a87, a91, a93, a94, a95, a96, a97, b1, b3, b4, b5, b6, b7, bp1, bp3, bp4, bp5, bp6, bp7, bp8, btilde1, btilde3, btilde4, btilde5, btilde6, btilde7, bptilde1, bptilde3, bptilde4, bptilde5, bptilde6, bptilde7, bptilde8, bptilde9 = cache.tab
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
    @inbounds ku[i] = uprev[i] + dt*(c5*duprev[i] + dt*(a61*k1.x[2][i] + a62*k2.x[2][i] + a63*k3.x[2][i] + a64*k4.x[2][i] + a65*k5.x[2][i]))
  end

  f.f2(t+dt*c5,ku,du,k6.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c6*duprev[i] + dt*(a71*k1.x[2][i] + a72*k2.x[2][i] + a73*k3.x[2][i] + a74*k4.x[2][i] + a75*k5.x[2][i] + a76*k6.x[2][i]))
  end

  f.f2(t+dt*c6,ku,du,k7.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c7*duprev[i] + dt*(a81*k1.x[2][i] + a82*k2.x[2][i] + a83*k3.x[2][i] + a84*k4.x[2][i] + a85*k5.x[2][i] + a86*k6.x[2][i] + a87*k7.x[2][i]))
  end

  f.f2(t+dt*c7,ku,du,k8.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c8*duprev[i] + dt*(a91*k1.x[2][i] + a93*k3.x[2][i] + a94*k4.x[2][i] + a95*k5.x[2][i] + a96*k6.x[2][i] + a97*k7.x[2][i])) # no a92 & a98
  end

  f.f2(t+dt*c8,ku,du,k9.x[2])
  @tight_loop_macros for i in uidx
    @inbounds u[i]  = uprev[i] + dt*(duprev[i] + dt*(b1 *k1.x[2][i] + b3 *k3.x[2][i] + b4 *k4.x[2][i] + b5 *k5.x[2][i] + b6 *k6.x[2][i] + b7 *k7.x[2][i])) # b1 -- b7, no b2
    @inbounds du[i] = duprev[i]+ dt*(bp1*k1.x[2][i] + bp3*k3.x[2][i] + bp4*k4.x[2][i] + bp5*k5.x[2][i] + bp6*k6.x[2][i] + bp7*k7.x[2][i] + bp8*k8.x[2][i]) # bp1 -- bp8, no bp2
  end

  f.f1(t+dt,u,du,k.x[1])
  f.f2(t+dt,u,du,k.x[2])
  if integrator.opts.adaptive
    uhat, duhat = utilde.x
    dtsq = dt^2
    @tight_loop_macros for i in uidx
      @inbounds uhat[i]  = dtsq*(btilde1*k1.x[2][i] + btilde3*k3.x[2][i] + btilde4*k4.x[2][i] + btilde5*k5.x[2][i] + btilde6*k6.x[2][i] + btilde7*k7.x[2][i])
      @inbounds duhat[i] = dt*(bptilde1*k1.x[2][i] + bptilde3*k3.x[2][i] + bptilde4*k4.x[2][i] + bptilde5*k5.x[2][i] + bptilde6*k6.x[2][i] + bptilde7*k7.x[2][i] + bptilde8*k8.x[2][i] + bptilde9*k9.x[2][i])
    end
    calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end

@muladd function perform_step!(integrator,cache::DPRKN12Cache,repeat_step=false)
  @unpack t,dt,f = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  @unpack tmp,atmp,fsalfirst,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k,utilde = cache
  @unpack c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, a21, a31, a32, a41, a42, a43, a51, a53, a54, a61, a63, a64, a65, a71, a73, a74, a75, a76, a81, a84, a85, a86, a87, a91, a93, a94, a95, a96, a97, a98, a101, a103, a104, a105, a106, a107, a108, a109, a111, a113, a114, a115, a116, a117, a118, a119, a1110, a121, a123, a124, a125, a126, a127, a128, a129, a1210, a1211, a131, a133, a134, a135, a136, a137, a138, a139, a1310, a1311, a1312, a141, a143, a144, a145, a146, a147, a148, a149, a1410, a1411, a1412, a1413, a151, a153, a154, a155, a156, a157, a158, a159, a1510, a1511, a1512, a1513, a1514, a161, a163, a164, a165, a166, a167, a168, a169, a1610, a1611, a1612, a1613, a1614, a1615, a171, a173, a174, a175, a176, a177, a178, a179, a1710, a1711, a1712, a1713, a1714, a1715, b1, b7, b8, b9, b10, b11, b12, b13, b14, b15, bp1, bp7, bp8, bp9, bp10, bp11, bp12, bp13, bp14, bp15, bp16, bp17, btilde1, btilde7, btilde8, btilde9, btilde10, btilde11, btilde12, btilde13, btilde14, btilde15, bptilde1, bptilde7, bptilde8, bptilde9, bptilde10, bptilde11, bptilde12, bptilde13, bptilde14, bptilde15, bptilde16, bptilde17 = cache.tab
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
    @inbounds ku[i] = uprev[i] + dt*(c4*duprev[i] + dt*(a51*k1.x[2][i] + a53*k3.x[2][i] + a54*k4.x[2][i])) # no a52
  end

  f.f2(t+dt*c4,ku,du,k5.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c5*duprev[i] + dt*(a61*k1.x[2][i] + a63*k3.x[2][i] + a64*k4.x[2][i] + a65*k5.x[2][i])) # no a62
  end

  f.f2(t+dt*c5,ku,du,k6.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c6*duprev[i] + dt*(a71*k1.x[2][i] + a73*k3.x[2][i] + a74*k4.x[2][i] + a75*k5.x[2][i] + a76*k6.x[2][i])) # no a72
  end

  f.f2(t+dt*c6,ku,du,k7.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c7*duprev[i] + dt*(a81*k1.x[2][i] + a84*k4.x[2][i] + a85*k5.x[2][i] + a86*k6.x[2][i] + a87*k7.x[2][i])) # no a82, a83
  end

  f.f2(t+dt*c7,ku,du,k8.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c8*duprev[i] + dt*(a91*k1.x[2][i] + a93*k3.x[2][i] + a94*k4.x[2][i] + a95*k5.x[2][i] + a96*k6.x[2][i] + a97*k7.x[2][i] + a98*k8.x[2][i])) # no a92
  end

  f.f2(t+dt*c8,ku,du,k9.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c9*duprev[i] + dt*(a101*k1.x[2][i] + a103*k3.x[2][i] + a104*k4.x[2][i] + a105*k5.x[2][i] + a106*k6.x[2][i] + a107*k7.x[2][i] + a108*k8.x[2][i] + a109*k9.x[2][i])) # no a102
  end

  f.f2(t+dt*c9,ku,du,k10.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c10*duprev[i] + dt*(a111*k1.x[2][i] + a113*k3.x[2][i] + a114*k4.x[2][i] + a115*k5.x[2][i] + a116*k6.x[2][i] + a117*k7.x[2][i] + a118*k8.x[2][i] + a119*k9.x[2][i] + a1110*k10.x[2][i])) # no a112
  end

  f.f2(t+dt*c10,ku,du,k11.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c11*duprev[i] + dt*(a121*k1.x[2][i] + a123*k3.x[2][i] + a124*k4.x[2][i] + a125*k5.x[2][i] + a126*k6.x[2][i] + a127*k7.x[2][i] + a128*k8.x[2][i] + a129*k9.x[2][i] + a1210*k10.x[2][i] + a1211*k11.x[2][i])) # no a122
  end

  f.f2(t+dt*c11,ku,du,k12.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c12*duprev[i] + dt*(a131*k1.x[2][i] + a133*k3.x[2][i] + a134*k4.x[2][i] + a135*k5.x[2][i] + a136*k6.x[2][i] + a137*k7.x[2][i] + a138*k8.x[2][i] + a139*k9.x[2][i] + a1310*k10.x[2][i] + a1311*k11.x[2][i] + a1312*k12.x[2][i])) # no a132
  end

  f.f2(t+dt*c12,ku,du,k13.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c13*duprev[i] + dt*(a141*k1.x[2][i] + a143*k3.x[2][i] + a144*k4.x[2][i] + a145*k5.x[2][i] + a146*k6.x[2][i] + a147*k7.x[2][i] + a148*k8.x[2][i] + a149*k9.x[2][i] + a1410*k10.x[2][i] + a1411*k11.x[2][i] + a1412*k12.x[2][i] + a1413*k13.x[2][i])) # no a142
  end

  f.f2(t+dt*c13,ku,du,k14.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c14*duprev[i] + dt*(a151*k1.x[2][i] + a153*k3.x[2][i] + a154*k4.x[2][i] + a155*k5.x[2][i] + a156*k6.x[2][i] + a157*k7.x[2][i] + a158*k8.x[2][i] + a159*k9.x[2][i] + a1510*k10.x[2][i] + a1511*k11.x[2][i] + a1512*k12.x[2][i] + a1513*k13.x[2][i] + a1514*k14.x[2][i])) # no a152
  end

  f.f2(t+dt*c14,ku,du,k15.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c15*duprev[i] + dt*(a161*k1.x[2][i] + a163*k3.x[2][i] + a164*k4.x[2][i] + a165*k5.x[2][i] + a166*k6.x[2][i] + a167*k7.x[2][i] + a168*k8.x[2][i] + a169*k9.x[2][i] + a1610*k10.x[2][i] + a1611*k11.x[2][i] + a1612*k12.x[2][i] + a1613*k13.x[2][i] + a1614*k14.x[2][i] + a1615*k15.x[2][i])) # no a162
  end

  f.f2(t+dt*c15,ku,du,k16.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = uprev[i] + dt*(c16*duprev[i] + dt*(a171*k1.x[2][i] + a173*k3.x[2][i] + a174*k4.x[2][i] + a175*k5.x[2][i] + a176*k6.x[2][i] + a177*k7.x[2][i] + a178*k8.x[2][i] + a179*k9.x[2][i] + a1710*k10.x[2][i] + a1711*k11.x[2][i] + a1712*k12.x[2][i] + a1713*k13.x[2][i] + a1714*k14.x[2][i] + a1715*k15.x[2][i])) # no a172, a1716
  end

  f.f2(t+dt*c16,ku,du,k17.x[2])
  @tight_loop_macros for i in uidx
    @inbounds u[i]  = uprev[i] + dt*(duprev[i] + dt*(b1 *k1.x[2][i] + b7 *k7.x[2][i]+ b8 *k8.x[2][i]+ b9 *k9.x[2][i]+ b10 *k10.x[2][i]+ b11 *k11.x[2][i]+ b12 *k12.x[2][i]+ b13 *k13.x[2][i]+ b14 *k14.x[2][i]+ b15 *k15.x[2][i])) # b1 & b7 -- b15
    @inbounds du[i] = duprev[i]+ dt*(bp1*k1.x[2][i] + bp7*k7.x[2][i]+ bp8*k8.x[2][i]+ bp9*k9.x[2][i]+ bp10*k10.x[2][i]+ bp11*k11.x[2][i]+ bp12*k12.x[2][i]+ bp13*k13.x[2][i]+ bp14*k14.x[2][i]+ bp15*k15.x[2][i]+ bp16*k16.x[2][i]+ bp17*k17.x[2][i]) # bp1 & bp7 -- bp17
  end

  f.f1(t+dt,u,du,k.x[1])
  f.f2(t+dt,u,du,k.x[2])
  if integrator.opts.adaptive
    uhat, duhat = utilde.x
    dtsq = dt^2
    @tight_loop_macros for i in uidx
      @inbounds uhat[i]  = dtsq*(btilde1*k1.x[2][i] + btilde7*k7.x[2][i] + btilde8*k8.x[2][i] + btilde9*k9.x[2][i] + btilde10*k10.x[2][i] + btilde11*k11.x[2][i] + btilde12*k12.x[2][i] + btilde13*k13.x[2][i] + btilde14*k14.x[2][i] + btilde15*k15.x[2][i]) # btilde1 & btilde7 -- btilde15
      @inbounds duhat[i] = dt*(bptilde1*k1.x[2][i] + bptilde7*k7.x[2][i] + bptilde8*k8.x[2][i] + bptilde9*k9.x[2][i] + bptilde10*k10.x[2][i] + bptilde11*k11.x[2][i] + bptilde12*k12.x[2][i] + bptilde13*k13.x[2][i] + bptilde14*k14.x[2][i] + bptilde15*k15.x[2][i] + bptilde16*k16.x[2][i] + bptilde17*k17.x[2][i]) # bptilde1 & bptilde7 -- bptilde17
    end
    calculate_residuals!(atmp, utilde, integrator.uprev, integrator.u, integrator.opts.abstol, integrator.opts.reltol)
    integrator.EEst = integrator.opts.internalnorm(atmp)
  end
end
