module OrdinaryDiffEqDefault

using OrdinaryDiffEq: Vern7, Vern8, Vern9, Vern6, Tsit5, Rosenbrock23, Rodas5P, FBDF, Krylov, KrylovFBDF

include("composite_algs.jl")
include("AutoSwitch.jl")
include("integrators/integrator_interface.jl")

end # module OrdinaryDiffEqDefault
