@doc raw"""
    FunctionMap(; scale_by_time = false)

A fixed timestep method for when the ODE is a discrete dynamical system. In the 
operator setting, this is equivalent to operator splitting for additive operators.

When `scale_by_time = true`, the method becomes `u_{n+1} = u_n + dt*f(u_n,p,t_n)`, 
otherwise it's `u_{n+1} = f(u_n,p,t_n)`.

!!! note
    This method requires a fixed timestep dt and is not adaptive.
"""
struct FunctionMap{scale_by_time} <: OrdinaryDiffEqAlgorithm end
FunctionMap(; scale_by_time = false) = FunctionMap{scale_by_time}()
