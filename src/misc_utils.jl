# Default nlsolve behavior, should move to FiniteDiff.jl

Base.@pure determine_chunksize(u,alg::DiffEqBase.DEAlgorithm) = determine_chunksize(u,get_chunksize(alg))
Base.@pure function determine_chunksize(u,CS)
  if CS != 0
    return CS
  else
    return ForwardDiff.pickchunksize(length(u))
  end
end

macro swap!(x,y)
  quote
    local tmp = $(esc(x))
    $(esc(x)) = $(esc(y))
    $(esc(y)) = tmp
  end
end

macro cache(expr)
  name = expr.args[2].args[1].args[1]
  fields = [x for x in expr.args[3].args if typeof(x)!=LineNumberNode]
  cache_vars = Expr[]
  jac_vars = Pair{Symbol,Expr}[]
  for x in fields
    if x.args[2] == :uType || x.args[2] == :rateType ||
       x.args[2] == :kType || x.args[2] == :uNoUnitsType
      push!(cache_vars,:(c.$(x.args[1])))
    elseif x.args[2] == :DiffCacheType
      push!(cache_vars,:(c.$(x.args[1]).du))
      push!(cache_vars,:(c.$(x.args[1]).dual_du))
    end
  end
  quote
    $expr
    $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
  end
end

# Nest one layer of value in order to get rid of possible Dual{Complex} or Complex{Dual} issues
# value should recurse for anything else.
function constvalue(::Type{T}) where T
  _T = DiffEqBase.value(T)
  return _T <: Complex ? DiffEqBase.value(real(_T)) : DiffEqBase.value(_T)
end
function constvalue(x)
  _x = DiffEqBase.value(x)
  return _x isa Complex ? DiffEqBase.value(real(_x)) : DiffEqBase.value(_x)
end

function diffdir(integrator::DiffEqBase.DEIntegrator)
  difference = maximum(abs, integrator.uprev)*sqrt(eps(typeof(integrator.t)))
  dir = integrator.tdir > zero(integrator.tdir) ?
          integrator.t > integrator.sol.prob.tspan[2] - difference ? -true :  true :
          integrator.t < integrator.sol.prob.tspan[2] + difference ?  true : -true
end

abstract type AbstractThreadingOption end
struct Sequential <: AbstractThreadingOption end
struct BaseThreads <: AbstractThreadingOption end
struct PolyesterThreads <: AbstractThreadingOption end

isthreaded(b::Bool) = b
isthreaded(::Sequential) = false
isthreaded(::BaseThreads) = true
isthreaded(::PolyesterThreads) = true

macro threaded(option, ex)
  quote
    opt = $(esc(option))
    if (opt === BaseThreads()) || ((opt isa Bool) && opt)
      $(esc(:(Threads.@threads $ex)))
    elseif opt === PolyesterThreads()
      $(esc(:(Polyester.@batch $ex)))
    else
      $(esc(ex))
    end
  end
end


