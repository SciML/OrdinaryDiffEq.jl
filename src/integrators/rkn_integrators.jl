## y'' = f(t, y, y')
## y(t₀) = y₀; y'(t₀) = y'₀
## kᵢ' = f(t₀+cᵢh, y₀+cᵢhy'₀+h²∑āᵢⱼk'ⱼ, y'₀+h∑aᵢⱼk'ⱼ)
## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
## y'₁ = y'₀ + h∑bᵢk'ᵢ

@inline function initialize!(integrator,cache::Nystrom4Cache,f=integrator.f)
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  uprev,duprev = integrator.uprev.x

  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Nystrom4Cache,f=integrator.f)
  @unpack t,dt,k = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  uidx = eachindex(integrator.uprev.x[1])
  @unpack tmp,fsalfirst,k₂,k₃,k₄,k = cache
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  k₁ = fsalfirst
  halfdt = dt/2
  dtsq = dt^2
  eighth_dtsq = dtsq/8
  half_dtsq = dtsq/2
  ttmp = t+halfdt

  f[2](ttmp,uprev,duprev,k₁.x[2])
  @tight_loop_macros for i in uidx
    ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
    @inbounds ku[i]  = @muladd uprev[i] + halfdt*duprev[i] + eighth_dtsq*k₁.x[2][i]
    ## y'₁ = y'₀ + h∑bᵢk'ᵢ
    @inbounds kdu[i] = @muladd duprev[i] + halfdt*k₁.x[2][i]
  end

  f[2](ttmp,ku,kdu,k₂.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i]  = @muladd uprev[i] + halfdt*duprev[i] + eighth_dtsq*k₁.x[2][i]
    @inbounds kdu[i] = @muladd duprev[i] + halfdt*k₂.x[2][i]
  end

  f[2](ttmp,ku,kdu,k₃.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i]  = @muladd uprev[i] + dt*duprev[i] + half_dtsq*k₃.x[2][i]
    @inbounds kdu[i] = @muladd duprev[i] + dt*k₃.x[2][i]
  end

  f[2](t+dt,ku,kdu,k₄.x[2])
  @tight_loop_macros for i in uidx
    @inbounds u[i]  = muladd(dt, duprev[i], muladd(dtsq/6, k₁.x[2][i] + k₂.x[2][i] + k₃.x[2][i], uprev[i]))
    @inbounds du[i] = muladd(dt/6, muladd(2, (k₂.x[2][i] + k₃.x[2][i]), k₁.x[2][i] + k₄.x[2][i]), duprev[i])
  end
  f[1](t+dt,ku,kdu,k.x[1])
  f[2](t+dt,ku,kdu,k.x[2])
end


@inline function initialize!(integrator,cache::Nystrom4VelocityIndependentCache,f=integrator.f)
  @unpack tmp,fsalfirst,k₂,k = cache
  uprev,duprev = integrator.uprev.x

  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Nystrom4VelocityIndependentCache,f=integrator.f)
  @unpack t,dt,k = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  uidx = eachindex(integrator.uprev.x[1])
  @unpack tmp,fsalfirst,k₂,k₃,k = cache
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  k₁ = fsalfirst
  halfdt = dt/2
  dtsq = dt^2
  eighth_dtsq = dtsq/8
  half_dtsq = dtsq/2
  ttmp = t+halfdt

  f[2](ttmp,uprev,duprev,k₁.x[2])
  @tight_loop_macros for i in uidx
    ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
    @inbounds ku[i] = @muladd uprev[i] + halfdt*duprev[i] + eighth_dtsq*k₁.x[2][i]
  end

  f[2](ttmp,ku,du,k₂.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = @muladd uprev[i] + dt*duprev[i] + half_dtsq*k₂.x[2][i]
  end

  f[2](t+dt,ku,du,k₃.x[2])
  @tight_loop_macros for i in uidx
    @inbounds u[i] = muladd(dt, duprev[i], muladd(dtsq/6, muladd(2, k₂.x[2][i], k₁.x[2][i]),uprev[i]))
    @inbounds du[i] = muladd(dt/6,muladd(4, k₂.x[2][i], k₁.x[2][i] + k₃.x[2][i]),duprev[i])
  end
  f[1](t+dt,ku,du,k.x[1])
  f[2](t+dt,ku,du,k.x[2])
end


@inline function initialize!(integrator,cache::Nystrom5VelocityIndependentCache,f=integrator.f)
  @unpack tmp,fsalfirst,k₂,k = cache
  uprev,duprev = integrator.uprev.x

  integrator.fsalfirst = fsalfirst
  integrator.fsallast = k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  f[1](integrator.t,uprev,duprev,integrator.k[2].x[1])
  f[2](integrator.t,uprev,duprev,integrator.k[2].x[2])
end

@inline function perform_step!(integrator,cache::Nystrom5VelocityIndependentCache,f=integrator.f)
  @unpack t,dt,k = integrator
  u,du = integrator.u.x
  uprev,duprev = integrator.uprev.x
  uidx = eachindex(integrator.uprev.x[1])
  @unpack tmp,fsalfirst,k₂,k₃,k = cache
  ku, kdu = integrator.cache.tmp.x[1], integrator.cache.tmp.x[2]
  k₁ = fsalfirst
  halfdt = dt/2
  dtsq = dt^2
  eighth_dtsq = dtsq/8
  half_dtsq = dtsq/2
  ttmp = t+halfdt

  f[2](ttmp,uprev,duprev,k₁.x[2])
  @tight_loop_macros for i in uidx
    ## y₁ = y₀ + hy'₀ + h²∑b̄ᵢk'ᵢ
    @inbounds ku[i] = @muladd uprev[i] + halfdt*duprev[i] + eighth_dtsq*k₁.x[2][i]
  end

  f[2](ttmp,ku,du,k₂.x[2])
  @tight_loop_macros for i in uidx
    @inbounds ku[i] = @muladd uprev[i] + dt*duprev[i] + half_dtsq*k₂.x[2][i]
  end

  f[2](t+dt,ku,du,k₃.x[2])
  @tight_loop_macros for i in uidx
    @inbounds u[i] = muladd(dt, duprev[i], muladd(dtsq/6, muladd(2, k₂.x[2][i], k₁.x[2][i]),uprev[i]))
    @inbounds du[i] = muladd(dt/6,muladd(4, k₂.x[2][i], k₁.x[2][i] + k₃.x[2][i]),duprev[i])
  end
  f[1](t+dt,ku,du,k.x[1])
  f[2](t+dt,ku,du,k.x[2])
end

