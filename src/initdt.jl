@muladd function ode_determine_initdt(u0,t,tdir,dtmax,abstol,reltol,internalnorm,prob::AbstractODEProblem{uType,tType,true},integrator) where {tType,uType}
  f = prob.f
  p = integrator.p
  oneunit_tType = oneunit(tType)
  dtmax_tdir = tdir*dtmax

  if eltype(u0) <: Number && !(typeof(integrator.alg) <: CompositeAlgorithm)
    cache = get_tmp_cache(integrator)
    sk = first(cache)
    @. sk = abstol+internalnorm(u0)*reltol
  else
    sk = @. abstol+internalnorm(u0)*reltol
  end

  if get_current_isfsal(integrator.alg, integrator.cache) && typeof(integrator) <: ODEIntegrator
    # Right now DelayDiffEq has issues with fsallast not being initialized
    f₀ = integrator.fsallast
    f(f₀,u0,p,t)
  else
    # TODO: use more caches
    f₀ = u0/t; fill!(f₀,zero(eltype(f₀)))
    f(f₀,u0,p,t)
  end

  # TODO: use more caches
  #tmp = cache[2]
  tmp = @. u0/sk

  d₀ = internalnorm(tmp)

  #=
  Try/catch around the linear solving. This will catch singular matrices defined
  by DAEs and thus we use the tType(1//10^(6)) default from Hairer. Note that
  this will not always catch singular matrices, an example from Andreas:

  julia> A = fill(rand(), 2, 2)
  2×2 Array{Float64,2}:
   0.637947  0.637947
   0.637947  0.637947

  julia> inv(A)
  2×2 Array{Float64,2}:
    9.0072e15  -9.0072e15
   -9.0072e15   9.0072e15

  The only way to make this more correct is to check

  issingular(A) = rank(A) < min(size(A)...)

  but that would introduce another svdfact in rank (which may not be possible
  anyways if the mass_matrix is not actually an array). Instead we stick to the
  user-chosen factorization. Sometimes this will cause `ftmp` to be absurdly
  large like shown there, but that later gets caught in the quick estimates
  below which then makes it spit out the default

  dt₀ = tType(1//10^(6))

  so that is a a cheaper way to get to the same place than an svdfact and
  still works for matrix-free definitions of the mass matrix.
  =#

  if prob.mass_matrix != I
    ftmp = similar(f₀)
    try
      integrator.alg.linsolve(ftmp, copy(prob.mass_matrix), f₀, true)
      f₀ .= ftmp
    catch
      return tType(1//10^(6))
    end
  end

  if any(isnan,f₀)
    warn("First function call produced NaNs. Exiting.")
  end

  @. tmp = f₀/sk*oneunit_tType
  d₁ = internalnorm(tmp)

  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  dt₀ = min(dt₀,dtmax_tdir)

  if typeof(one(tType)) <: AbstractFloat && dt₀ < 10eps(tType)*oneunit(tType)
    # This catches Andreas' non-singular example
    # should act like it's singular
    return tdir*tType(1//10^(6))
  end

  dt₀_tdir = tdir*dt₀

  u₁ = similar(u0) # required by DEDataArray
  @. u₁ = u0 + dt₀_tdir*f₀
  f₁ = similar(f₀)
  f(f₁,u₁,p,t+dt₀_tdir)

  if prob.mass_matrix != I
    integrator.alg.linsolve(ftmp, prob.mass_matrix, f₁, false)
    f₁ .= ftmp
  end

  @. tmp = (f₁-f₀)/sk*oneunit_tType
  d₂ = internalnorm(tmp)/dt₀*oneunit_tType
  # Hairer has d₂ = sqrt(sum(abs2,tmp))/dt₀, note the lack of norm correction

  max_d₁d₂ = max(d₁,d₂)
  if max_d₁d₂ <= 1//Int64(10)^(15)
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max_d₁d₂))/get_current_alg_order(integrator.alg,integrator.cache)))
  end
  dt = tdir*min(100dt₀,dt₁,dtmax_tdir)
end

@muladd function ode_determine_initdt(u0,t,tdir,dtmax,abstol,reltol,internalnorm,prob::AbstractODEProblem{uType,tType,false},integrator) where {uType,tType}
  f = prob.f
  p = prob.p
  oneunit_tType = oneunit(tType)
  dtmax_tdir = tdir*dtmax

  sk = abstol+internalnorm(u0)*reltol
  d₀ = internalnorm(u0/sk)

  f₀ = f(u0,p,t)
  if any((isnan(x) for x in f₀))
    error("First function call produced NaNs. Exiting.")
  end

  d₁ = internalnorm(f₀/sk*oneunit_tType)

  if d₀ < 1//10^(5) || d₁ < 1//10^(5)
    dt₀ = tType(1//10^(6))
  else
    dt₀ = tType((d₀/d₁)/100)
  end
  dt₀ = min(dt₀,dtmax_tdir)
  dt₀_tdir = tdir*dt₀

  u₁ = u0 + dt₀_tdir*f₀
  f₁ = f(u₁,p,t+dt₀_tdir)

  d₂ = internalnorm((f₁-f₀)./sk*oneunit_tType)/dt₀*oneunit_tType

  max_d₁d₂ = max(d₁, d₂)
  if max_d₁d₂ <= 1//Int64(10)^(15)
    dt₁ = max(tType(1//10^(6)),dt₀*1//10^(3))
  else
    dt₁ = tType(10.0^(-(2+log10(max_d₁d₂))/get_current_alg_order(integrator.alg,integrator.cache)))
  end
  dt = tdir*min(100dt₀,dt₁,dtmax_tdir)
end
