type RHS_IIF1_Scalar{F,CType,tType} <: Function
  f::F
  t::tType
  dt::tType
  tmp::CType
end

function (p::RHS_IIF1_Scalar)(u,resid)
  resid[1] = u[1] - p.tmp - p.dt*p.f[2](p.t+p.dt,u[1])[1]
end

type RHS_IIF2_Scalar{F,CType,tType} <: Function
  f::F
  t::tType
  dt::tType
  tmp::CType
end

function (p::RHS_IIF2_Scalar)(u,resid)
  resid[1] = u[1] - p.tmp - 0.5p.dt*p.f[2](p.t+p.dt,u[1])[1]
end

function initialize!(integrator,cache::Union{IIF1ConstantCache,IIF2ConstantCache})
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  A = integrator.f[1](integrator.t,integrator.u)
  cache.uhold[1] = f[2](integrator.t,integrator.uprev)
  integrator.fsalfirst = A*integrator.uprev .+ cache.uhold[1]

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsalfirst
end

function perform_step!(integrator,cache::Union{IIF1ConstantCache,IIF2ConstantCache},repeat_step=false)
  @unpack t,dt,uprev,u,f = integrator
  @unpack uhold,rhs,nl_rhs = cache

  # If adaptive, this should be computed after and cached
  A = integrator.f[1]
  if typeof(cache) <: IIF1ConstantCache
    tmp = expm(A*dt)*(uprev)
  elseif typeof(cache) <: IIF2ConstantCache
    @muladd tmp = expm(A*dt)*(@. uprev + 0.5dt*uhold[1]) # This uhold only works for non-adaptive
  end

  if integrator.success_iter > 0 && !integrator.u_modified
    uhold[1] = current_extrapolant(t+dt,integrator)
  end # else uhold is previous value.

  rhs.t = t
  rhs.dt = dt
  rhs.tmp = tmp
  nlres = integrator.alg.nlsolve(nl_rhs,uhold)
  uhold[1] = integrator.f[2](t+dt,nlres[1])
  u = nlres[1]
  integrator.fsallast = A*u .+ uhold[1]
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

type RHS_IIF1{F,uType,tType,DiffCacheType,SizeType,uidxType} <: Function
  f::F
  tmp::uType
  t::tType
  dt::tType
  dual_cache::DiffCacheType
  sizeu::SizeType
  uidx::uidxType
end
function (p::RHS_IIF1)(u,resid)
  du = get_du(p.dual_cache, eltype(u))
  p.f[2](p.t+p.dt,reshape(u,p.sizeu),du)
  @. resid = u - p.tmp - p.dt*du
end

type RHS_IIF2{F,uType,tType,DiffCacheType,SizeType,uidxType} <: Function
  f::F
  tmp::uType
  t::tType
  dt::tType
  dual_cache::DiffCacheType
  sizeu::SizeType
  uidx::uidxType
end
function (p::RHS_IIF2)(u,resid)
  du = get_du(p.dual_cache, eltype(u))
  p.f[2](p.t+p.dt,reshape(u,p.sizeu),du)
  @. resid = u - p.tmp - 0.5p.dt*du
end

function initialize!(integrator,cache::Union{IIF1Cache,IIF2Cache})
  integrator.fsalfirst = cache.fsalfirst
  integrator.fsallast = cache.k
  integrator.kshortsize = 2
  integrator.k = eltype(integrator.sol.k)(integrator.kshortsize)
  A = integrator.f[1]
  integrator.f[2](integrator.t,integrator.uprev,cache.rtmp1)
  A_mul_B!(cache.k,A,integrator.uprev)
  @. integrator.fsalfirst = cache.k + cache.rtmp1
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

function perform_step!(integrator,cache::Union{IIF1Cache,IIF2Cache},repeat_step=false)
  @unpack rtmp1,tmp,k = cache
  @unpack uhold,rhs,nl_rhs = cache
  @unpack t,dt,uprev,u,f = integrator

  @. k = uprev
  if typeof(cache) <: IIF2Cache
    @muladd @. k = k + 0.5dt*rtmp1
  end

  A = integrator.f[1]
  M = expm(A*dt)
  A_mul_B!(tmp,M,k)

  if integrator.success_iter > 0 && !integrator.u_modified
    current_extrapolant!(uhold,t+dt,integrator)
  end # else uhold is previous value.

  rhs.t = t
  rhs.dt = dt
  rhs.tmp = tmp
  rhs.uidx = eachindex(u)
  rhs.sizeu = size(u)
  nlres = integrator.alg.nlsolve(nl_rhs,uhold)

  copy!(u,nlres)
  integrator.f[2](t+dt,nlres,rtmp1)
  integrator.fsallast .= A*u .+ rtmp1
end
