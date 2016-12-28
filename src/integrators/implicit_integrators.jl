function ode_solve{uType<:Number,tType,ksEltype,F,rateType,O}(integrator::ODEIntegrator{ImplicitEuler,uType,tType,ksEltype,F,rateType,O})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}
  function rhs_ie(u,resid,u_old,t,dt)
    resid[1] = u[1] - u_old[1] - dt*f(t+dt,u)[1]
  end
  uhold::Vector{uType} = Vector{uType}(1)
  u_old::Vector{uType} = Vector{uType}(1)
  uhold[1] = u; u_old[1] = u
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      u_old[1] = uhold[1]
      nlres = NLsolve.nlsolve((uhold,resid)->rhs_ie(uhold,resid,u_old,t,dt),uhold,autodiff=autodiff)
      uhold[1] = nlres.zero[1]
      if integrator.opts.calck
        k = f(t+dt,uhold[1])
      end
      u = uhold[1]
      @ode_loopfooter
    end
  end
  @ode_postamble
end

function ode_solve{uType<:AbstractArray,tType,ksEltype,F,rateType,O}(integrator::ODEIntegrator{ImplicitEuler,uType,tType,ksEltype,F,rateType,O})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}
  uidx = eachindex(u)
  if autodiff
    dual_cache = DiffCache(u)
    rhs_ie = (u,resid,u_old,t,dt,dual_cache) -> begin
      du = get_du(dual_cache, eltype(u))
      f(t+dt,reshape(u,sizeu),du)
      for i in uidx
        resid[i] = u[i] - u_old[i] - dt*du[i]
      end
    end
  else
    dual_cache = similar(u)
    rhs_ie = (u,resid,u_old,t,dt,du) -> begin
      f(t+dt,reshape(u,sizeu),du)
      for i in uidx
        resid[i] = u[i] - u_old[i] - dt*du[i]
      end
    end
  end

  uhold = vec(u); u_old = similar(u)
  cache = (u,u_old,dual_cache,uprev,kprev)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      copy!(u_old,uhold)
      nlres = NLsolve.nlsolve((uhold,resid)->rhs_ie(uhold,resid,u_old,t,dt,dual_cache),uhold,autodiff=autodiff)
      uhold[:] = nlres.zero
      if integrator.opts.calck
        f(t+dt,u,k)
      end
      @ode_loopfooter
    end
  end
  u = reshape(uhold,sizeu...)
  @ode_postamble
end

function ode_solve{uType<:AbstractArray,tType,ksEltype,F,rateType,O}(integrator::ODEIntegrator{Trapezoid,uType,tType,ksEltype,F,rateType,O})
  @ode_preamble
  local nlres::NLsolve.SolverResults{uEltype}
  uidx = eachindex(u)
  if autodiff
    cache1 = DiffCache(u)
    cache2 = DiffCache(u)
    dto2 = dt/2
    rhs_trap = (u,resid,u_old,t,dt,cache1,cache2) -> begin
      du1 = get_du(cache1, eltype(u)); du2 = get_du(cache2, eltype(u_old))
      f(t,reshape(u_old,sizeu),du2)
      f(t+dt,reshape(u,sizeu),du1)
      for i in uidx
        resid[i] = u[i] - u_old[i] - dto2*(du1[i]+du2[i])
      end
    end
  else
    cache1 = similar(u)
    cache2 = similar(u)
    dto2 = dt/2
    rhs_trap = (u,resid,u_old,t,dt,du1,du2) -> begin
      f(t,reshape(u_old,sizeu),du2)
      f(t+dt,reshape(u,sizeu),du1)
      for i in uidx
        resid[i] = u[i] - u_old[i] - dto2*(du1[i]+du2[i])
      end
    end
  end
  uhold = vec(u); u_old = similar(u)

  cache = (u,u_old,cache1,cache2,uprev,kprev)
  @inbounds for T in Ts
    while t < T
      @ode_loopheader
      copy!(u_old,uhold)
      nlres = NLsolve.nlsolve((uhold,resid)->rhs_trap(uhold,resid,u_old,t,dt,cache1,cache2),uhold,autodiff=autodiff)
      uhold[:] = nlres.zero
      if integrator.opts.calck
        f(t+dt,u,k)
      end
      @ode_loopfooter
    end
  end
  u = reshape(uhold,sizeu...)
  @ode_postamble
end

function ode_solve{uType<:Number,tType,ksEltype,F,rateType,O}(integrator::ODEIntegrator{Trapezoid,uType,tType,ksEltype,F,rateType,O})
  @ode_preamble
  dto2::tType = dt/2
  function rhs_trap(u,resid,u_old,t,dt)
    resid[1] = u[1] - u_old[1] - dto2*(f(t,u_old)[1] + f(t+dt,u)[1])
  end
  local nlres::NLsolve.SolverResults{uEltype}
  uhold::Vector{uType} = Vector{uType}(1)
  u_old::Vector{uType} = Vector{uType}(1)
  uhold[1] = u; u_old[1] = u
  @inbounds for T in Ts
      while t < T
      @ode_loopheader
      u_old[1] = uhold[1]
      nlres = NLsolve.nlsolve((uhold,resid)->rhs_trap(uhold,resid,u_old,t,dt),uhold,autodiff=autodiff)
      uhold[1] = nlres.zero[1]
      if integrator.opts.calck
        k = f(t+dt,uhold[1])
      end
      u = uhold[1]
      @ode_loopfooter
    end
  end
  @ode_postamble
end
