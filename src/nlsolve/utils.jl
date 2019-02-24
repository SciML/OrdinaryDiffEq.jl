"""
    qrdelete!(Q, R, k)

Delete the left-most column of F = Q[:, 1:k] * R[1:k, 1:k] by updating Q and R.
Only Q[:, 1:(k-1)] and R[1:(k-1), 1:(k-1)] are valid on exit.
"""
function qrdelete!(Q::AbstractMatrix, R::AbstractMatrix, k::Int)
  n, m = size(Q)
  m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
  1 ≤ k ≤ m || throw(ArgumentError())

  # apply Givens rotations
  for i in 2:k
      g = first(givens(R, i - 1, i, i))
      lmul!(g, R)
      rmul!(Q, g')
  end

  # move columns of R
  @inbounds for j in 1:(k-1)
    for i in 1:(k-1)
      R[i, j] = R[i, j + 1]
    end
  end

  Q, R
end

"""
    qradd!(Q, R, v, k)

Replace the right-most column of F = Q[:, 1:k] * R[1:k, 1:k] with v by updating Q and R.
This implementation modifies vector v as well. Only Q[:, 1:k] and R[1:k, 1:k] are valid on
exit.
"""
function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::AbstractVector, k::Int)
  n, m = size(Q)
  n == length(v) || throw(DimensionMismatch())
  m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
  1 ≤ k ≤ m || throw(ArgumentError())

  @inbounds for i in 1:(k-1)
    q = view(Q, :, i)
    r = dot(q, v)

    R[i, k] = r
    axpy!(-r, q, v)
  end

  @inbounds begin
    d = norm(v)
    R[k, k] = d
    @. Q[:, k] = v / d
  end

  Q, R
end

set_new_W!(nlsolver::NLSolver, val::Bool) = set_new_W!(nlsolver.cache, val)
set_new_W!(nlcache::NLNewtonCache, val::Bool) = (nlcache.new_W = val; nothing)
set_new_W!(nlcache::AbstractNLSolverCache, val::Bool) = nothing

isnewton(nlsolver::NLSolver) = isnewton(nlsolver.cache)
isnewton(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = true
isnewton(nlcache::AbstractNLSolverCache) = false

get_W(nlsolver::NLSolver) = get_W(nlsolver.cache)
get_W(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = nlcache.W

function get_κtol(nlalg::Union{NLAnderson,NLFunctional,NLNewton}, uTolType, reltol)
  κ = nlalg.κ === nothing ? uTolType(1//100) : uTolType(nlalg.κ)
  tol = nlalg.tol === nothing ? uTolType(min(0.03, first(reltol)^(0.5))) : uTolType(nlalg.tol)
  κ * tol
end

DiffEqBase.@def iipnlsolve begin
  # define additional fields of cache of non-linear solver
  z = similar(u); dz = similar(u); tmp = similar(u); b = similar(u)
  k = zero(rate_prototype)

  # adapt options of non-linear solver to current integration problem
  uTolType = real(uBottomEltypeNoUnits)
  κtol = get_κtol(alg.nlsolve, uTolType, reltol)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)

    # check if `nf` is linear
    islin = f isa Union{ODEFunction,SplitFunction} && islinear(nf.f)

    if islin
      # get the operator
      J = nf.f
      W = WOperator(f.mass_matrix, dt, J, true)
    else
      if DiffEqBase.has_jac(f) && !DiffEqBase.has_invW(f) && f.jac_prototype !== nothing
        W = WOperator(f, dt, true)
        J = nothing # is J = W.J better?
      else
        J = fill(zero(uEltypeNoUnits),length(u),length(u)) # uEltype?
        W = similar(J)
      end
    end

    nlcache = NLNewtonCache(true,W)
  elseif alg.nlsolve isa NLFunctional
    z₊ = similar(z)

    nlcache = NLFunctionalCache(z₊)
  elseif alg.nlsolve isa NLAnderson
    z₊ = similar(z)

    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = [zero(z) for i in 1:max_history]
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)
    dzold = zero(z)
    z₊old = zero(z)

    nlcache = NLAndersonCache(z₊,dzold,z₊old,Δz₊s,Q,R,γs,alg.nlsolve.aa_start,alg.nlsolve.droptol)
  end

  # create non-linear solver
  nlsolver = NLSolver{true,typeof(z),typeof(k),uTolType,typeof(γ),typeof(c),typeof(nlcache)}(z,dz,tmp,b,k,one(uTolType),κtol,γ,c,alg.nlsolve.max_iter,10000,nlcache)

  # define additional fields of cache
  fsalfirst = zero(rate_prototype)
  if alg.nlsolve isa NLNewton
    if islin
      du1 = rate_prototype
      uf = nothing
      jac_config = nothing
      linsolve = alg.linsolve(Val{:init},nf,u)
    else
      du1 = zero(rate_prototype)
      # if the algorithm specializes on split problems the use `nf`
      uf = DiffEqDiffTools.UJacobianWrapper(nf,t,p)
      jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,tmp,dz)
      linsolve = alg.linsolve(Val{:init},uf,u)
    end
  else
    J = nothing
    W = nothing
    du1 = rate_prototype
    uf = nothing
    jac_config = nothing
    linsolve = nothing
  end
end

DiffEqBase.@def oopnlsolve begin
  # define additional fields of cache of non-linear solver (all aliased)
  z = uprev; dz = z; tmp = z; b = z; k = rate_prototype

  # define tolerances
  uTolType = real(uBottomEltypeNoUnits)
  κtol = get_κtol(alg.nlsolve, uTolType, reltol)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)

    # only use `nf` if the algorithm specializes on split eqs
    uf = DiffEqDiffTools.UDerivativeWrapper(nf,t,p)

    islin = f isa Union{ODEFunction,SplitFunction} && islinear(nf.f)
    if islin || DiffEqBase.has_jac(f)
      # get the operator
      J = islin ? nf.f : f.jac(uprev, p, t)
      if !isa(J, DiffEqBase.AbstractDiffEqLinearOperator)
        J = DiffEqArrayOperator(J)
      end
      W = WOperator(f.mass_matrix, dt, J, false)
    else
      if DiffEqBase.has_jac(f)
        J = f.jac(uprev, p, t)
      else
        if alg_autodiff(alg)
          J = jacobian_autodiff(uf, uprev)
        else
          J = jacobian_finitediff(uf, uprev, alg.diff_type)
        end
      end
      W = J isa Number ? J : lu(J; check=false)
    end

    nlcache = NLNewtonConstantCache(W)
  elseif alg.nlsolve isa NLFunctional
    uf = nothing

    nlcache = NLFunctionalConstantCache()
  elseif alg.nlsolve isa NLAnderson
    uf = nothing

    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = Vector{typeof(z)}(undef, max_history)
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    nlcache = NLAndersonConstantCache(Δz₊s,Q,R,γs,alg.nlsolve.aa_start,alg.nlsolve.droptol)
  end

  # create non-linear solver
  nlsolver = NLSolver{false,typeof(z),typeof(k),uTolType,typeof(γ),typeof(c),typeof(nlcache)}(z,dz,tmp,b,k,one(uTolType),κtol,γ,c,alg.nlsolve.max_iter,10000,nlcache)
end
