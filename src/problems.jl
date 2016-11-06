"""
`ODEProblem`

Wraps the data which defines an ODE problem

```math
\\frac{du}{dt} = f(t,u)
```

with initial condition ``u₀``.

### Constructors

`ODEProblem(f,u₀;analytic=nothing)` : Defines the ODE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the ODE.
* `u₀`: The initial condition.
* `analytic`: A function which describes the solution.
* `knownanalytic`: True if the solution is given.
* `numvars`: The number of variables in the system

"""
type ODEProblem{uType,uEltype,tType} <: AbstractODEProblem
  f::Function
  u₀::uType
  analytic::Function
  knownanalytic::Bool
  numvars::Int
  isinplace::Bool
  tspan::Vector{tType}
end

function ODEProblem(f::Function,u₀,tspan=[0,1.];analytic=nothing)
  isinplace = numparameters(f)>=3
  if analytic==nothing
    knownanalytic = false
    analytic=(t,u,du)->0
  else
    knownanalytic = true
  end
  if typeof(u₀) <: Number
    sizeu = (1,)
    numvars = 1
  else
    numvars = length(u₀)[end]
  end
  ODEProblem{typeof(u₀),eltype(u₀)}(f,u₀,analytic,knownanalytic,numvars,isinplace,tspan)
end
