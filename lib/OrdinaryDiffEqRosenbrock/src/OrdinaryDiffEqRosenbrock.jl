module OrdinaryDiffEqRosenbrock

import OrdinaryDiffEq: OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
                       OrdinaryDiffEqRosenbrockAlgorithm, rosenbrock_docstring,
                       rosenbrock_wanner_docstring,
                       OrdinaryDiffEqMutableCache, alg_order, alg_adaptive_order,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm
using DiffEqBase, MuladdMacro

include("algorithms.jl")
include("alg_utils.jl")
include("rosenbrock_tableaus.jl")
include("generic_rosenbrock.jl")
include("rosenbrock_caches.jl")
include("integrator_interface.jl")
include("integrator_utils.jl")
include("interp_func.jl")
include("stiff_addsteps.jl")
include("rosenbrock_interpolants.jl")
include("rosenbrock_tableaus.jl")

export Rosenbrock23, Rosenbrock32, RosShamp4, Veldd4, Velds4, GRK4T, GRK4A,
       Ros4LStab, ROS3P, Rodas3, Rodas23W, Rodas3P, Rodas4, Rodas42, Rodas4P, Rodas4P2,
       Rodas5, Rodas5P, Rodas5Pe, Rodas5Pr,
       RosenbrockW6S4OS, ROS34PW1a, ROS34PW1b, ROS34PW2, ROS34PW3, ROS34PRw, ROS3PRL,
       ROS3PRL2, ROK4a,
       ROS2, ROS2PR, ROS2S, ROS3, ROS3PR, Scholz4_7

end
