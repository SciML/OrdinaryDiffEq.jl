@muladd function ode_determine_initdt{tType,uType}(u0,t::tType,tdir,dtmax,abstol,reltol,internalnorm,prob::AbstractODEProblem{uType,tType,true},order)
  f = prob.f
  f₀ = zeros(u0./t); f₁ = zeros(u0./t); u₁ = zeros(u0); sk = zeros(u0);
  # Hack to  make a generic u0 with no units, https://github.com/JuliaLang/julia/issues/22216
  typeof(u0[1]) <: AbstractArray ? tmp = zeros(u0,typeof(ones(u0[1]))) : tmp = zeros(u0,typeof(one(u0[1])))

  @. sk = abstol+abs(u0)*reltol
  @. tmp = u0/sk
  d₀ = internalnorm(tmp)

  f(t,u0,f₀)
  if any((isnan(x) for x in f₀))
    warn("First function call produced NaNs. Exiting.")
  end

  @. tmp = f₀/sk*oneunit(tType)
  d₁ = internalnorm(tmp)

  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  dt₀ = min(dt₀,tdir*dtmax)

  @. u₁ = u0 + tdir*dt₀*f₀
  f(t+tdir*dt₀,u₁,f₁)

  @. tmp = (f₁-f₀)/(abstol+abs(u0)*reltol)*oneunit(tType)
  d₂ = internalnorm(tmp)/dt₀*oneunit(tType)
  # Hairer has d₂ = sqrt(sum(abs2,tmp))/dt₀, note the lack of norm correction

  if max(d₁,d₂) <= 1//Int64(10)^(15)
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)))/order))
  end
  dt = tdir*min(100dt₀,dt₁,tdir*dtmax)
end

@muladd function ode_determine_initdt{uType,tType}(u0::uType,t,tdir,dtmax,abstol,reltol,internalnorm,prob::AbstractODEProblem{uType,tType,false},order)
  f = prob.f

  sk = @. abstol+abs(u0)*reltol
  d₀ = internalnorm(@. u0/sk*oneunit(tType))

  f₀ = f(t,u0)
  if any((isnan(x) for x in f₀))
    error("First function call produced NaNs. Exiting.")
  end

  d₁ = internalnorm(@. f₀/sk*oneunit(tType))

  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  dt₀ = min(dt₀,tdir*dtmax)

  u₁ = @. u0 + tdir*dt₀*f₀
  f₁ = f(t+tdir*dt₀,u₁)

  d₂ = internalnorm(@. (f₁-f₀)/(abstol+abs(u0)*reltol)*oneunit(tType))/dt₀*oneunit(tType)

  if max(d₁,d₂) <= 1//Int64(10)^(15)
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max(d₁,d₂)))/order))
  end
  dt = tdir*min(100dt₀,dt₁,tdir*dtmax)
end
