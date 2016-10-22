module OrdinaryDiffEq

using DiffEqBase
import DiffEqBase: solve
using IterativeSolvers, Parameters, GenericSVD, ForwardDiff, InplaceOps,
      Sundials, ParameterizedFunctions, Ranges, NLsolve

import Base: linspace
import Plots: plot
import ForwardDiff.Dual

macro def(name, definition)
    return quote
        macro $name()
            esc($(Expr(:quote, definition)))
        end
    end
end


typealias KW Dict{Symbol,Any}
AbstractArrayOrVoid = Union{AbstractArray,Void}
NumberOrVoid = Union{Number,Void}
FunctionOrVoid = Union{Function,Void}

#Constants

const TEST_FLOPS_CUTOFF = 1e10
const initialized_backends = Set{Symbol}()

include("backends.jl")
include("sensitivity.jl")
include("misc_utils.jl")

include("problems.jl")
include("solutions.jl")

include("solve/ode_integrators.jl")
include("solve/ode_constants.jl")
include("solve/ode_callbacks.jl")
include("solve/ode_solve.jl")
include("solve/ode_dense.jl")
include("solve/unrolled_tableaus.jl")

include("premade_problems.jl")

#Types
export ODESolution, ODEProblem


#ODE Example Problems
export prob_ode_linear, prob_ode_bigfloatlinear, prob_ode_2Dlinear,
       prob_ode_large2Dlinear, prob_ode_bigfloat2Dlinear, prob_ode_rigidbody,
       prob_ode_2Dlinear_notinplace, prob_ode_vanderpol, prob_ode_vanderpol_stiff,
       prob_ode_lorenz, prob_ode_rober, prob_ode_threebody

 #General Functions
 export solve

 #Callback Necessary
 export ode_addsteps!, ode_interpolant, DIFFERENTIALEQUATIONSJL_SPECIALDENSEALGS,
        @ode_callback, @ode_event, @ode_change_cachesize, @ode_change_deleteat,
        @ode_terminate, @ode_savevalues, copyat_or_push!

 export constructDP5, constructVern6, constructDP8, constructDormandPrince, constructFeagin10,
        constructFeagin12, constructFeagin14


end # module
