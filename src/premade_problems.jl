using OrdinaryDiffEq, ParameterizedFunctions
srand(100)

### ODE Examples

# Linear ODE
f = (t,u) -> (1.01*u)
analytic = (t,u₀) -> u₀*exp(1.01*t)
"""
Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u₀=1/2``, ``α=1.01``, and solution

```math
u(t) = u₀e^{αt}
```

with Float64s
"""
prob_ode_linear = ODEProblem(f,1/2,analytic=analytic)

const linear_bigα = parse(BigFloat,"1.01")
f = (t,u) -> (linear_bigα*u)
analytic = (t,u₀) -> u₀*exp(linear_bigα*t)
"""
Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u₀=1/2``, ``α=1.01``, and solution

```math
u(t) = u₀e^{αt}
```

with BigFloats
"""
prob_ode_bigfloatlinear = ODEProblem(f,parse(BigFloat,"0.5"),analytic=analytic)

f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = 1.01*u[i]
  end
end
analytic = (t,u₀) -> u₀*exp.(1.01*t)
"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u₀=1/2``, ``α=1.01``, and solution

```math
u(t) = u₀e^{αt}
```

with Float64s
"""
prob_ode_2Dlinear = ODEProblem(f,rand(4,2),analytic=analytic)
"""
100x100 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u₀=1/2``, ``α=1.01``, and solution

```math
u(t) = u₀e^{αt}
```

with Float64s
"""
prob_ode_large2Dlinear = ODEProblem(f,rand(100,100),analytic=analytic)

f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u₀=1/2``, ``α=1.01``, and solution

```math
u(t) = u₀e^{αt}
```

with BigFloats
"""
prob_ode_bigfloat2Dlinear = ODEProblem(f,map(BigFloat,rand(4,2)).*ones(4,2)/2,analytic=analytic)
f = (t,u) -> 1.01*u
"""
4x2 version of the Linear ODE

```math
\\frac{du}{dt} = αu
```

with initial condition ``u₀=1/2``, ``α=1.01``, and solution

```math
u(t) = u₀e^{αt}
```

on Float64. Purposefully not in-place as a test.
"""
prob_ode_2Dlinear_notinplace = ODEProblem(f,rand(4,2),analytic=analytic)

## Lotka-Volterra


f = @ode_def LotkaVolterra begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=1.5 b=1 c=3 d=1

"""
Lotka-Voltera Equations

```math
\\frac{dx}{dt} = ax - bxy
\\frac{dy}{dt} = -cy + dxy
```

with initial condition ``x=y=1``
"""
prb_ode_lotkavoltera = ODEProblem(f,[1;1])

## Fitzhugh-Nagumo

f = @ode_def FitzhughNagumo begin
  dv = v - v^3/3 -w + l
  dw = τinv*(v +  a - b*w)
end a=0.7 b=0.8 τinv=(1/12.5) l=0.5
"""
Fitzhugh-Nagumo

```math
\\frac{dv}{dt} = v - \\frac{v^3}{3} - w + I_{est}
τ \\frac{dw}{dt} = v + a -bw
```

with initial condition ``v=w=1``
"""
prob_ode_fitzhughnagumo = ODEProblem(f,[1;1])

#Van der Pol Equations
f = @ode_def VanDerPol begin
  dy = μ*(1-x^2)*y - x
  dx = 1*y
end μ=>1.

"""
Van der Pol Equations

```math
\\begin{align}
\\frac{dx}{dt} &= y \\\\
\\frac{dy}{dt} &= μ(1-x^2)y -x
\\end{align}
```

with ``μ=1.0`` and ``u₀=[0,\\sqrt{3}]``

Non-stiff parameters.
"""
prob_ode_vanderpol = ODEProblem(f,[0;sqrt(3)])

f = VanDerPol(μ=1e6)
"""Van der Pol Equations

```math
\\begin{align}
\\frac{dx}{dt} &= y \\\\
\\frac{dy}{dt} &= μ(1-x^2)y -x
\\end{align}
```

with ``μ=10^6`` and ``u₀=[0,\\sqrt{3}]``

Stiff parameters.
"""
prob_ode_vanderpol_stiff = ODEProblem(f,[0;sqrt(3)])

# ROBER

f = @ode_def Rober begin
  dy₁ = -k₁*y₁+k₃*y₂*y₃
  dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
  dy₃ =  k₂*y₂^2
end k₁=>0.04 k₂=>3e7 k₃=>1e4

"""
The Robertson biochemical reactions:

```math
\\begin{align}
\\frac{dy₁}{dt} &= -k₁y₁+k₃y₂y₃  \\\\
\\frac{dy₂}{dt} &=  k₁y₁-k₂y₂^2-k₃y₂y₃ \\\\
\\frac{dy₃}{dt} &=  k₂y₂^2
\\end{align}
```

where ``k₁=0.04``, ``k₂=3\\times10^7``, ``k₃=10^4``. For details, see:

Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 129

Usually solved on `[0,1e11]`
"""
prob_ode_rober = ODEProblem(f,[1.0;0.0;0.0])

# Three Body
const threebody_μ = parse(BigFloat,"0.012277471"); const threebody_μ′ = 1 - threebody_μ

f = (t,u,du) -> begin
  # 1 = y₁
  # 2 = y₂
  # 3 = y₁'
  # 4 = y₂'
  D₁ = ((u[1]+threebody_μ)^2 + u[2]^2)^(3/2)
  D₂ = ((u[1]-threebody_μ′)^2 + u[2]^2)^(3/2)
  du[1] = u[3]
  du[2] = u[4]
  du[3] = u[1] + 2u[4] - threebody_μ′*(u[1]+threebody_μ)/D₁ - threebody_μ*(u[1]-threebody_μ′)/D₂
  du[4] = u[2] - 2u[3] - threebody_μ′*u[2]/D₁ - threebody_μ*u[2]/D₂
end
"""
The ThreeBody problem as written by Hairer:

```math
\\begin{align}
y₁′′ &= y₁ + 2y₂′ - μ′\\frac{y₁+μ}{D₁} - μ\\frac{y₁-μ′}{D₂} \\\\
y₂′′ &= y₂ - 2y₁′ - μ′\\frac{y₂}{D₁} - μ\\frac{y₂}{D₂} \\\\
D₁ &= ((y₁+μ)^2 + y₂^2)^{3/2} \\\\
D₂ &= ((y₁-μ′)^2+y₂^2)^{3/2} \\\\
μ &= 0.012277471 \\\\
μ′ &=1-μ
\\end{align}
```

From Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 129

Usually solved on `t₀ = 0.0`; `T = parse(BigFloat,"17.0652165601579625588917206249")`
Periodic with that setup.
"""
prob_ode_threebody = ODEProblem(f,[0.994, 0.0, 0.0, parse(BigFloat,"-2.00158510637908252240537862224")])

# Rigid Body Equations

f = @ode_def RigidBody begin
  dy₁  = I₁*y₂*y₃
  dy₂  = I₂*y₁*y₃
  dy₃  = I₃*y₁*y₂
end I₁=>-2 I₂=>1.25 I₃=>-.5

"""
Rigid Body Equations

```math
\\begin{align}
\\frac{dy₁}{dt}  &= I₁y₂y₃ \\\\
\\frac{dy₂}{dt}  &= I₂y₁y₃ \\\\
\\frac{dy₃}{dt}  &= I₃y₁y₂
\\end{align}
```

with ``I₁=-2``, ``I₂=1.25``, and ``I₃=-1/2``.

The initial condition is ``y=[1.0;0.0;0.9]``.

From Solving Differential Equations in R by Karline Soetaert

or Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 244

Usually solved from 0 to 20.
"""
prob_ode_rigidbody = ODEProblem(f,[1.0,0.0,0.9])

# Pleiades Problem

f = (t,u,du) -> begin
  x = view(u,1:7)   # x
  y = view(u,8:14)  # y
  v = view(u,15:21) # x′
  w = view(u,22:28) # y′
  du[1:7] .= v
  du[8:14].= w
  for i in 14:21
    du[i] = zero(u)
  end
  for i=1:7,j=1:7
    if i != j
      r = ((x[i]-x[j])^2 + (y[i] - y[j])^2)^(3/2)
      du[14+i] += j*(x[j] - x[i])/r
      du[21+i] += j*(y[j] - y[i])/r
    end
  end
end
"""
Pleides Problem

```math
\\begin{align}
xᵢ′′ &= \\sum_{j≠i} mⱼ(xⱼ-xᵢ)/rᵢⱼ \\\\
yᵢ′′ &= \\sum_{j≠i} mⱼ(yⱼ-yᵢ)/rᵢⱼ
\\end{align}
```

where

```math
rᵢⱼ = ((xᵢ-xⱼ)^2 + (yᵢ-yⱼ)^2)^{3/2}
```

and inital condtions are

```math
\\begin{align}
x₁(0)&=3  \\\\
x₂(0)&=3  \\\\
x₃(0)&=-1  \\\\
x₄(0)&=-3  \\\\
x₅(0)&=2  \\\\
x₆(0)&=-2  \\\\
x₇(0)&=2  \\\\
y₁(0)&=3  \\\\
y₂(0)&=-3  \\\\
y₃(0)&=2  \\\\
y₄(0)&=0  \\\\
y₅(0)&=0  \\\\
y₆(0)&=-4  \\\\
y₇(0)&=4
\\end{align}
```

and with ``xᵢ′(0)=yᵢ′(0)=0`` except for

```math
\\begin{align}
x₆′(0)&=1.75 \\\\
x₇′(0)&=-1.5 \\\\
y₄′(0)&=-1.25 \\\\
y₅′(0)&=1
\\end{align}
```

From Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 244

Usually solved from 0 to 3.
"""
prob_ode_pleides = ODEProblem(f,[3.0,3.0,-1.0,-3.0,2.0,-2.0,2.0,3.0,-3.0,2.0,0,0,-4.0,4.0,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
