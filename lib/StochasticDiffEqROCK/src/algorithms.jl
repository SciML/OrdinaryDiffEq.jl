struct SROCK1{interpretation, E} <: StochasticDiffEqAlgorithm
    eigen_est::E
end
function SROCK1(; interpretation = SciMLBase.AlgorithmInterpretation.Ito, eigen_est = nothing)
    return SROCK1{interpretation, typeof(eigen_est)}(eigen_est)
end

# Weak Order 2

for Alg in [:SROCK2, :KomBurSROCK2, :SROCKC2]
    @eval begin
        struct $Alg{E} <: StochasticDiffEqAlgorithm
            eigen_est::E
        end
        $Alg(; eigen_est = nothing) = $Alg(eigen_est)
    end
end

@doc """
    SROCK2(;eigen_est=nothing)

**SROCK2: Second-Order Stabilized Runge-Kutta Chebyshev Method**

Second-order stabilized explicit method with weak order 2.0 for mildly stiff SDE problems.

## Method Properties
- **Strong Order**: 1.0
- **Weak Order**: 2.0
- **Time stepping**: Fixed step size with extended stability
- **Stability**: Extended along negative real axis

## When to Use
- When higher accuracy than SROCK1 is needed
- Parabolic PDEs requiring better weak convergence
- Problems where second-order accuracy justifies increased cost

## References
- Second-order ROCK methods for stochastic problems
""" SROCK2

@doc """
    KomBurSROCK2(;eigen_est=nothing)

**KomBurSROCK2: Komori-Burrage Second-Order SROCK Method**

Alternative second-order stabilized method with different coefficients and stability properties.

## Method Properties
- **Strong Order**: 1.0
- **Weak Order**: 2.0
- **Time stepping**: Fixed step size with extended stability
- **Stability**: Extended along negative real axis

## When to Use
- Alternative to SROCK2 with different stability characteristics
- When SROCK2 performance is unsatisfactory
- Benchmarking against other second-order ROCK methods

## References
- Komori and Burrage stabilized methods
""" KomBurSROCK2

@doc """
    SROCKC2(;eigen_est=nothing)

**SROCKC2: Conservative Second-Order SROCK Method**

Conservative second-order stabilized method designed for robust performance.

## Method Properties
- **Strong Order**: 1.0
- **Weak Order**: 2.0
- **Time stepping**: Fixed step size with extended stability
- **Stability**: Conservative stability region, more robust

## When to Use
- When robustness is more important than efficiency
- For difficult problems where other ROCK methods fail
- As a fallback option for problematic cases

## References
- Conservative ROCK methods for stochastic problems
""" SROCKC2

# ROCK stabilization for EM
"""
    SROCKEM(;strong_order_1=true, eigen_est=nothing)

**SROCKEM: ROCK-Stabilized Euler-Maruyama Method**

Fixed step Euler-Maruyama method with first-order ROCK stabilization for handling stiff problems.

## Method Properties

  - **Strong Order**: 1.0 (default) or 0.5 (if `strong_order_1=false`)
  - **Weak Order**: 1.0 (default) or 0.5 (if `strong_order_1=false`)
  - **Time stepping**: Fixed step size with ROCK stabilization
  - **Noise types**: 1-dimensional, diagonal, and multi-dimensional noise
  - **SDE interpretation**: Ito only
  - **Stability**: ROCK stabilization for moderate stiffness

## Parameters

  - `strong_order_1::Bool = true`: Use strong/weak order 1.0 (true) or 0.5 (false)
  - `eigen_est`: Eigenvalue estimation for stability (automatic if `nothing`)

## When to Use

  - Stiff problems where standard EM fails
  - When ROCK stabilization is preferred over full implicit treatment
  - Problems requiring Euler-Maruyama structure with enhanced stability
  - Multi-dimensional stiff SDEs

## Algorithm Description

Combines Euler-Maruyama discretization with ROCK stabilization techniques to extend the stability region without requiring linear solves.

## References

  - ROCK stabilization techniques applied to SDEs
  - Stabilized Euler methods for stiff problems
"""
struct SROCKEM{E} <: StochasticDiffEqAlgorithm
    strong_order_1::Bool
    eigen_est::E
end
SROCKEM(; strong_order_1 = true, eigen_est = nothing) = SROCKEM(strong_order_1, eigen_est)
"""
    SKSROCK(;post_processing=false, eigen_est=nothing)

**SKSROCK: SK-SROCK Stabilized Method**

Fixed step stabilized explicit method for stiff Ito problems with enhanced stability domain and optional post-processing.

## Method Properties

  - **Strong Order**: 0.5 (up to 2.0 with post-processing)
  - **Weak Order**: 1.0 (up to 2.0 with post-processing)
  - **Time stepping**: Fixed step size with enhanced stability
  - **Noise types**: 1-dimensional, diagonal, and multi-dimensional noise
  - **SDE interpretation**: Ito only
  - **Stability**: Better stability domain than SROCK1

## Parameters

  - `post_processing::Bool = false`: Enable post-processing for higher accuracy (experimental)
  - `eigen_est`: Eigenvalue estimation for stability (automatic if `nothing`)

## When to Use

  - Stiff Ito problems requiring better stability than SROCK1
  - Ergodic dynamical systems (with post-processing)
  - Problems where enhanced stability domain is crucial
  - When experimenting with post-processing techniques

## Post-Processing (Experimental)

  - Can achieve order 2 accuracy for ergodic systems
  - Particularly useful for Brownian dynamics
  - Currently under development - use with caution

## Algorithm Features

  - Enhanced stability compared to SROCK1
  - Handles various noise structures
  - Optional post-processing for specialized applications

## References

  - SK-SROCK methods for stochastic problems
  - Post-processing techniques for ergodic systems
"""
struct SKSROCK{E} <: StochasticDiffEqAlgorithm
    post_processing::Bool
    eigen_est::E
end
function SKSROCK(; post_processing = false, eigen_est = nothing)
    return SKSROCK(post_processing, eigen_est)
end
"""
    TangXiaoSROCK2(;version_num=5, eigen_est=nothing)

**TangXiaoSROCK2: Tang-Xiao Second-Order SROCK Method**

Fixed step size stabilized explicit method with multiple variants offering different stability domains.

## Method Properties

  - **Strong Order**: 1.0
  - **Weak Order**: 2.0
  - **Time stepping**: Fixed step size with extended stability
  - **Noise types**: Various (depends on version)
  - **SDE interpretation**: Ito only
  - **Stability**: Version-dependent stability domains

## Parameters

  - `version_num::Int = 5`: Choose version 1-5 with different stability characteristics
  - `eigen_est`: Eigenvalue estimation for stability (automatic if `nothing`)

## When to Use

  - When experimenting with different stability domains
  - Problems requiring weak order 2.0 with fixed steps
  - Benchmarking different ROCK variants
  - **Note**: Currently under development

## Versions

  - Versions 1-5 offer different stability domains
  - Version 5 (default) typically provides good general performance
  - Choose version based on problem-specific stability requirements

## Development Status

  - Method is under active development
  - Use with caution in production code
  - Consider more established ROCK methods for critical applications

## References

  - Tang and Xiao, "Second-order SROCK methods for stochastic problems"
"""
struct TangXiaoSROCK2{E} <: StochasticDiffEqAlgorithm
    version_num::Int
    eigen_est::E
end
function TangXiaoSROCK2(; version_num = 5, eigen_est = nothing)
    return TangXiaoSROCK2(version_num, eigen_est)
end
