struct DefaultInit <: DiffEqBase.DAEInitializationAlgorithm end
struct NoInit <: DiffEqBase.DAEInitializationAlgorithm end

struct ShampineCollocationInit{T} <: DiffEqBase.DAEInitializationAlgorithm
  initdt::T
end
ShampineCollocationInit() = ShampineCollocationInit(nothing)

struct BrownFullBasicInit{T} <: DiffEqBase.DAEInitializationAlgorithm
  abstol::T
end
BrownFullBasicInit() = BrownFullBasicInit(1e-10)

## Notes

#=
differential_vars = [any(!iszero,x) for x in eachcol(M)]

A column should be zero for an algebraic variable, since that means that the
derivative term doesn't show up in any equations (i.e. is an algebraic variable).
The rows are not necessarily non-zero, for example a flux condition between two
differential variables. But if it's a condition that doesn't involve the algebraic
variable, then the system is not Index 1!

=#

## Expansion

function DiffEqBase.initialize_dae!(integrator::ODEIntegrator,
                  initializealg = integrator.initializealg)
  _initialize_dae!(integrator, integrator.sol.prob,
           initializealg,
           Val(DiffEqBase.isinplace(integrator.sol.prob)))
end

## Default algorithms

function _initialize_dae!(integrator, prob::ODEProblem,
             alg::DefaultInit, x::Val{true})
  _initialize_dae!(integrator, prob,
          BrownFullBasicInit(), x)
end

function _initialize_dae!(integrator, prob::ODEProblem,
             alg::DefaultInit, x::Val{false})
  _initialize_dae!(integrator, prob,
          BrownFullBasicInit(), x)
end

function _initialize_dae!(integrator, prob::DAEProblem,
             alg::DefaultInit, x::Val{false})
  if prob.differential_vars === nothing
    _initialize_dae!(integrator, prob,
            ShampineCollocationInit(), x)
  else
    _initialize_dae!(integrator, prob,
            BrownFullBasicInit(), x)
  end
end

function _initialize_dae!(integrator, prob::DAEProblem,
             alg::DefaultInit, x::Val{true})
  if prob.differential_vars === nothing
    _initialize_dae!(integrator, prob,
            ShampineCollocationInit(), x)
  else
    _initialize_dae!(integrator, prob,
            BrownFullBasicInit(), x)
  end
end

## NoInit

function _initialize_dae!(integrator, prob::ODEProblem,
             alg::NoInit, x::Val{true})
end

function _initialize_dae!(integrator, prob::ODEProblem,
             alg::NoInit, x::Val{false})
end

function _initialize_dae!(integrator, prob::DAEProblem,
             alg::NoInit, x::Val{false})
end

function _initialize_dae!(integrator, prob::DAEProblem,
             alg::NoInit, x::Val{true})
end

## ShampineCollocationInit

#=
The method:

du = (u-u0)/h
Solve for `u`

=#

function _initialize_dae!(integrator, prob::ODEProblem, alg::ShampineCollocationInit, ::Val{true})
  @unpack p, t, f = integrator
  M = integrator.f.mass_matrix
  dtmax = integrator.opts.dtmax
  tmp = first(get_tmp_cache(integrator))
  u0 = integrator.u

  dt = if alg.initdt === nothing
    integrator.dt != 0 ? min(integrator.dt/5, dtmax) : 1//1000 # Haven't implemented norm reduction
  else
    alg.initdt
  end

  algebraic_vars = [all(iszero,x) for x in eachcol(M)]
  algebraic_eqs  = [all(iszero,x) for x in eachrow(M)]
  (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
  update_coefficients!(M, u0, p, t)
  f(tmp, u0, p, t)
  tmp .= ArrayInterface.restructure(tmp,algebraic_eqs .* vec(tmp))

  integrator.opts.internalnorm(tmp, t) <= integrator.opts.abstol && return

  if isdefined(integrator.cache, :nlsolver)
    # backward Euler
    nlsolver = integrator.cache.nlsolver
    oldγ, oldc, oldmethod, olddt = nlsolver.γ, nlsolver.c, nlsolver.method, integrator.dt
    nlsolver.tmp .= integrator.uprev
    nlsolver.γ, nlsolver.c = 1, 1
    nlsolver.method = DIRK
    integrator.dt = dt
    z = nlsolve!(nlsolver, integrator, integrator.cache)
    nlsolver.γ, nlsolver.c, nlsolver.method, integrator.dt = oldγ, oldc, oldmethod, olddt
    # TODO: failure handling
    nlsolvefail(nlsolver) && @warn "ShampineCollocationInit DAE initialization algorithm failed with dt=$dt. Try to adjust initdt like `ShampineCollocationInit(initdt)`."
    @.. integrator.u = integrator.uprev + z
  else
    isad = alg_autodiff(integrator.alg)
    _tmp = isad ? DiffEqBase.dualcache(tmp, Val{ForwardDiff.pickchunksize(length(tmp))}) : tmp
    nlequation! = @closure (out,u) -> begin
      update_coefficients!(M,u,p,t)
      #M * (u-u0)/dt - f(u,p,t)
      tmp = isad ? DiffEqBase.get_tmp(_tmp, u) : _tmp
      @. tmp = (u - u0)/dt
      mul!(vec(out),M,vec(tmp))
      f(tmp,u,p,t)
      out .-= tmp
      nothing
    end
    r = nlsolve(nlequation!, u0, autodiff=isad ? :forward : :central, method = :newton)
    integrator.u .= r.zero
  end
  recursivecopy!(integrator.uprev,integrator.u)
  if alg_extrapolates(integrator.alg)
    recursivecopy!(integrator.uprev2,integrator.uprev)
  end
  return nothing
end

function _initialize_dae!(integrator, prob::ODEProblem, alg::ShampineCollocationInit, ::Val{false})
  @unpack p, t, f = integrator
  u0 = integrator.u
  M = integrator.f.mass_matrix
  dtmax = integrator.opts.dtmax

  dt = if alg.initdt === nothing
    integrator.dt != 0 ? min(integrator.dt/5, dtmax) : 1//1000 # Haven't implemented norm reduction
  else
    alg.initdt
  end

  algebraic_vars = [all(iszero,x) for x in eachcol(M)]
  algebraic_eqs  = [all(iszero,x) for x in eachrow(M)]
  (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
  update_coefficients!(M,u0,p,t)
  du = f(u0,p,t)
  resid = @view vec(du)[algebraic_eqs]

  integrator.opts.internalnorm(resid,t) <= integrator.opts.abstol && return

  if isdefined(integrator.cache, :nlsolver)
    # backward Euler
    nlsolver = integrator.cache.nlsolver
    oldγ, oldc, oldmethod, olddt = nlsolver.γ, nlsolver.c, nlsolver.method, integrator.dt
    nlsolver.tmp .= integrator.uprev
    nlsolver.γ, nlsolver.c = 1, 1
    nlsolver.method = DIRK
    integrator.dt = dt
    z = nlsolve!(nlsolver, integrator, integrator.cache)
    nlsolver.γ, nlsolver.c, nlsolver.method, integrator.dt = oldγ, oldc, oldmethod, olddt
    # TODO: failure handling
    nlsolvefail(nlsolver) && @warn "ShampineCollocationInit DAE initialization algorithm failed with dt=$dt. Try to adjust initdt like `ShampineCollocationInit(initdt)`."
    @.. integrator.u = integrator.uprev + z
  else
    nlequation_oop = @closure u -> begin
      update_coefficients!(M,u,p,t)
      M * (u-u0)/dt - f(u,p,t)
    end

    nlequation! = @closure (out,u) -> out .= nlequation_oop(u)

    integrator.u = nlsolve(nlequation!, u0).zero
  end
  integrator.uprev = integrator.u
  if alg_extrapolates(integrator.alg)
    integrator.uprev2 = integrator.uprev
  end
  return
end

function _initialize_dae!(integrator, prob::DAEProblem,
             alg::ShampineCollocationInit, ::Val{true})

  @unpack p, t, f = integrator
  u0 = integrator.u

  dtmax = integrator.opts.dtmax
  tmp = get_tmp_cache(integrator)[1]
  resid = get_tmp_cache(integrator)[2]

  dt = t != 0 ? min(t/1000,dtmax) : dtmax # Haven't implemented norm reduction

   nlequation! = @closure (out,u) -> begin
     #M * (u-u0)/dt - f(u,p,t)
    @. tmp = (u - u0)/dt
    f(out,tmp,u,p,t)
    nothing
   end

  nlequation!(tmp,u0)
  f(resid,integrator.du,u0,p,t)
   integrator.opts.internalnorm(resid,t) <= integrator.opts.abstol && return

  integrator.u .= nlsolve(nlequation!, u0).zero
  recursivecopy!(integrator.uprev,integrator.u)
  if alg_extrapolates(integrator.alg)
    recursivecopy!(integrator.uprev2,integrator.uprev)
  end
  return
end

function _initialize_dae!(integrator, prob::DAEProblem,
             alg::ShampineCollocationInit, ::Val{false})

  @unpack p, t, f = integrator
  u0 = integrator.u
  dtmax = integrator.opts.dtmax

  dt = t != 0 ? min(t/1000,dtmax/10) : dtmax # Haven't implemented norm reduction

  nlequation_oop = u -> begin
    f((u-u0)/dt,u,p,t)
  end

  nlequation! = @closure (out,u) -> out .= nlequation_oop(u)

  resid = f(integrator.du,u0,p,t)
  integrator.opts.internalnorm(resid,t) <= integrator.opts.abstol && return

  integrator.u = nlsolve(nlequation!, u0).zero
  integrator.uprev = integrator.u
  if alg_extrapolates(integrator.alg)
    integrator.uprev2 = integrator.uprev
  end
end

## BrownFullBasic

#=
The method:

Keep differential variables constant
Solve for the algebraic variables

=#

function _initialize_dae!(integrator, prob::ODEProblem, alg::BrownFullBasicInit, ::Val{true})
  @unpack p, t, f = integrator
  u = integrator.u
  M = integrator.f.mass_matrix
  update_coefficients!(M,u,p,t)
  algebraic_vars = [all(iszero,x) for x in eachcol(M)]
  algebraic_eqs  = [all(iszero,x) for x in eachrow(M)]
  (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return
  tmp = get_tmp_cache(integrator)[1]

  f(tmp,u,p,t)

  tmp .= ArrayInterface.restructure(tmp,algebraic_eqs .* vec(tmp))

  integrator.opts.internalnorm(tmp,t) <= alg.abstol && return
  alg_u = @view u[algebraic_vars]

  isad = alg_autodiff(integrator.alg)
  if isad
    chunk = Val{ForwardDiff.pickchunksize(count(algebraic_vars))}
    _tmp = DiffEqBase.dualcache(tmp, chunk)
    _du_tmp = DiffEqBase.dualcache(tmp, chunk)
  else
    _tmp, _du_tmp = tmp, similar(tmp)
  end
  nlequation = @closure (out, x) -> begin
    uu = isad ? DiffEqBase.get_tmp(_tmp, x) : _tmp
    du_tmp = isad ? DiffEqBase.get_tmp(_du_tmp, x) : _du_tmp
    copyto!(uu, integrator.u)
    alg_uu = @view uu[algebraic_vars]
    alg_uu .= x
    f(du_tmp, uu, p, t)
    out .= @view du_tmp[algebraic_eqs]
    return nothing
  end

  r = nlsolve(nlequation, u[algebraic_vars], autodiff=isad ? :forward : :central, method=:newton)
  alg_u .= r.zero

  recursivecopy!(integrator.uprev,integrator.u)
  if alg_extrapolates(integrator.alg)
    recursivecopy!(integrator.uprev2,integrator.uprev)
  end

  return
end

function _initialize_dae!(integrator, prob::ODEProblem,
             alg::BrownFullBasicInit, ::Val{false})
  @unpack p, t, f = integrator

  u0 = integrator.u
  M = integrator.f.mass_matrix
  update_coefficients!(M,u0,p,t)
  algebraic_vars = [all(iszero,x) for x in eachcol(M)]
  algebraic_eqs  = [all(iszero,x) for x in eachrow(M)]
  (iszero(algebraic_vars) || iszero(algebraic_eqs)) && return

  du = f(u0,p,t)
  resid = @view vec(du)[algebraic_eqs]

  integrator.opts.internalnorm(resid,t) <= alg.abstol && return

  if u0 isa Number
    # This doesn't fix static arrays!
    u = [u0]
  else
    u = u0
  end

  alg_u = @view u[algebraic_vars]

  nlequation = @closure (out,x) -> begin
    alg_u .= x
    du = f(u,p,t)
    out .= @view du[algebraic_eqs]
  end
  r = nlsolve(nlequation, u0[algebraic_vars])
  alg_u .= r.zero

  if u0 isa Number
    # This doesn't fix static arrays!
    integrator.u = first(u)
  else
    integrator.u = u
  end

  integrator.uprev = integrator.u
  if alg_extrapolates(integrator.alg)
    integrator.uprev2 = integrator.uprev
  end

  return
end

function _initialize_dae!(integrator, prob::DAEProblem,
             alg::BrownFullBasicInit, ::Val{true})
  @unpack p, t, f = integrator
  differential_vars = prob.differential_vars
  u = integrator.u
  du = integrator.du

  tmp = get_tmp_cache(integrator)[1]
  f(tmp, du, u, p, t)

  if integrator.opts.internalnorm(tmp,t) <= alg.abstol
    return
  elseif differential_vars === nothing
    error("differential_vars must be set for DAE initialization to occur. Either set consistent initial conditions, differential_vars, or use a different initialization algorithm.")
  end

  nlequation = @closure (out, x) -> begin
    @. du = ifelse(differential_vars,x,du)
    @. u  = ifelse(differential_vars,u,x)
    f(out, du, u, p, t)
  end

  r = nlsolve(nlequation, ifelse.(differential_vars,du,u))

  @. du = ifelse(differential_vars,r.zero,du)
  @. u  = ifelse(differential_vars,u,r.zero)

  recursivecopy!(integrator.uprev,integrator.u)
  if alg_extrapolates(integrator.alg)
    recursivecopy!(integrator.uprev2,integrator.uprev)
  end

  return
end

function _initialize_dae!(integrator, prob::DAEProblem,
              alg::BrownFullBasicInit, ::Val{false})
  @unpack p, t, f = integrator
  differential_vars = prob.differential_vars

  if integrator.opts.internalnorm(f(integrator.du, integrator.u, p, t),t) <= alg.abstol
    return
  elseif differential_vars === nothing
    error("differential_vars must be set for DAE initialization to occur. Either set consistent initial conditions, differential_vars, or use a different initialization algorithm.")
  end

  if integrator.u isa Number && integrator.du isa Number
    # This doesn't fix static arrays!
    u = [integrator.u]
    du = [integrator.du]
  else
    u = integrator.u
    du = integrator.du
  end

  nlequation = @closure (out,x) -> begin
    @. du = ifelse(differential_vars,x,du)
    @. u  = ifelse(differential_vars,u,x)
    out .= f(du, u, p, t)
  end

  r = nlsolve(nlequation, ifelse.(differential_vars,du,u))

  @. du = ifelse(differential_vars,r.zero,du)
  @. u  = ifelse(differential_vars,u,r.zero)

  if integrator.u isa Number && integrator.du isa Number
    # This doesn't fix static arrays!
    integrator.u = first(u)
    integrator.du = first(du)
  else
    integrator.u = u
    integrator.du = du
  end

  integrator.uprev = integrator.u
  if alg_extrapolates(integrator.alg)
    integrator.uprev2 = integrator.uprev
  end

  return
end
