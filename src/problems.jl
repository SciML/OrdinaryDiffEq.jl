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
* `numvars`: The number of variables in the system

"""
type ODEProblem{uType,tType} <: AbstractODEProblem
  f::Base.Callable
  u0::uType
  isinplace::Bool
  tspan::Vector{tType}
end

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
* `numvars`: The number of variables in the system

"""
type ODETestProblem{uType,tType} <: AbstractODEProblem
  f::Base.Callable
  u0::uType
  analytic::Base.Callable
  isinplace::Bool
  tspan::Vector{tType}
end

function ODEProblem(f::Base.Callable,u0,tspan)
  isinplace = numparameters(f)>=3
  ODEProblem{typeof(u0),eltype(tspan)}(f,u0,isinplace,tspan)
end

function ODETestProblem(f::Base.Callable,u0,analytic,tspan=[0,1])
  isinplace = numparameters(f)>=3
  ODETestProblem{typeof(u0),eltype(tspan)}(f,u0,analytic,isinplace,tspan)
end
