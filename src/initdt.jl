function ode_determine_initdt{tType,uType}(u0,t::tType,tdir,dtmax,abstol,reltol,internalnorm,prob::AbstractODEProblem{uType,tType,true},order)
  f = prob.f
  f₀ = zeros(u0./t); f₁ = zeros(u0./t); u₁ = zeros(u0); sk = zeros(u0); tmp = zeros(u0)
  uidx = eachindex(u0)
  #sk = abstol+abs.(u0).*reltol
  #d₀ = internalnorm(u0./sk)

  @tight_loop_macros for i in uidx
    @inbounds sk[i] = (abstol+abs(u0[i])*reltol)
    @inbounds tmp[i] = u0[i]/sk[i]
  end
  d₀ = internalnorm(tmp)

  f(t,u0,f₀)
  #d₁ = internalnorm((f₀./sk*tType(1))/tType(1))

  @tight_loop_macros for i in uidx
    tmp[i] = (f₀[i]./sk[i]*tType(1))/tType(1)
  end
  d₁ = internalnorm(tmp)

  T0 = typeof(d₀)
  T1 = typeof(d₁)
  if d₀ < T0(1//10^(5)) || d₁ < T1(1//10^(5))
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  dt₀ = min(dt₀,tdir*dtmax)
  
  #@. u₁ = @muladd u0 + tdir*dt₀*f₀
  @tight_loop_macros for i in uidx
    @inbounds u₁[i] = u0[i] + tdir*dt₀*f₀[i]
  end

  f(t+tdir*dt₀,u₁,f₁)

  #tmp = (f₁.-f₀)./(abstol+abs.(u0).*reltol)*tType(1)
  @tight_loop_macros for i in uidx
    tmp[i] = (f₁[i]-f₀[i])/(abstol+abs(u0[i])*reltol)*tType(1)
  end

  d₂ = internalnorm(tmp)/dt₀
  # Hairer has d₂ = sqrt(sum(abs2,tmp))/dt₀, note the lack of norm correction
  unitless_max = max(d₁/typeof(d₁)(one(d₁)),d₂/typeof(d₂)(one(d₂)))
  if unitless_max<=1//10^(15)
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(unitless_max))/(order)))
  end
  dt = tdir*min(100dt₀,dt₁,tdir*dtmax)
end

function ode_determine_initdt{uType,tType}(u0::uType,t,tdir,dtmax,abstol,reltol,internalnorm,prob::AbstractODEProblem{uType,tType,false},order)
  f = prob.f
  sk = abstol+abs.(u0).*reltol
  d₀ = internalnorm(u0./sk)
  f₀ = f(t,u0)
  d₁ = internalnorm(f₀./sk)
  T0 = typeof(d₀)
  T1 = typeof(d₁)
  if d₀ < T0(1//10^(5)) || d₁ < T1(1//10^(5))
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  dt₀ = min(dt₀,tdir*dtmax)
  u₁ = u0 + tdir*dt₀*f₀
  f₁ = f(t+tdir*dt₀,u₁)
  d₂ = internalnorm((f₁-f₀)./(abstol+abs.(u0).*reltol))/dt₀*tType(1)
  if max(d₁,d₂) <= T1(1//10^(15))
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)/T1(1)))/(order)))
  end
  dt = tdir*min(100dt₀,dt₁,tdir*dtmax)
end
