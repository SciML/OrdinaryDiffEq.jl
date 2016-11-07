"""
`ODEProblem`

Wraps the data which defines an ODE problem

```math
\\frac{du}{dt} = f(t,u)
```

with initial condition ``u0``.

### Constructors

`ODEProblem(f,u0;analytic=nothing)` : Defines the ODE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the ODE.
* `u0`: The initial condition.
* `analytic`: A function which describes the solution.
* `knownanalytic`: True if the solution is given.
* `numvars`: The number of variables in the system

"""
type ODEProblem{uType,uEltype,tType} <: AbstractODEProblem
  f::Base.Callable
  u0::uType
  analytic::Base.Callable
  knownanalytic::Bool
  numvars::Int
  isinplace::Bool
  tspan::Vector{tType}
end

function ODEProblem(f::Function,u0,tspan=[0,1];analytic=nothing)
  isinplace = numparameters(f)>=3
  if analytic==nothing
    knownanalytic = false
    analytic=(t,u,du)->0
  else
    knownanalytic = true
  end
  if typeof(u0) <: Number
    sizeu = (1,)
    numvars = 1
  else
    numvars = length(u0)[end]
  end
  ODEProblem{typeof(u0),eltype(u0),eltype(tspan)}(f,u0,analytic,knownanalytic,numvars,isinplace,tspan)
end
