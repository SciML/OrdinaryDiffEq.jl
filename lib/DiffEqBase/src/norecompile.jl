struct OrdinaryDiffEqTag end

const NORECOMPILE_ARGUMENT_MESSAGE = """
No-recompile mode is only supported for state arguments
of type `Vector{Float64}`, time arguments of `Float64`
and parameter arguments of type `Vector{Float64}` or
`SciMLBase.NullParameters`.
"""

struct NoRecompileArgumentError <: Exception
    args::Any
end

function Base.showerror(io::IO, e::NoRecompileArgumentError)
    println(io, NORECOMPILE_ARGUMENT_MESSAGE)
    print(io, "Attempted arguments: ")
    return print(io, e.args)
end

function unwrap_fw(fw::FunctionWrapper)
    return fw.obj[]
end

# Default dispatch assumes no ForwardDiff, gets added in the new dispatch
function wrapfun_iip(ff, inputs)
    return FunctionWrappersWrappers.FunctionWrappersWrapper(
        Void(ff), (typeof(inputs),), (Nothing,)
    )
end

# 3-arg fallback: when ForwardDiff extension is not loaded, ignore chunk size
wrapfun_iip(ff, inputs, ::Val) = wrapfun_iip(ff, inputs)

function wrapfun_oop(ff, inputs)
    return FunctionWrappersWrappers.FunctionWrappersWrapper(
        ff, (typeof(inputs),), (typeof(inputs[1]),)
    )
end

# Wrap an in-place Jacobian function jac!(J, u, p, t) -> Nothing.
# Unlike the RHS, the Jacobian is not called with Dual numbers
# (the analytical Jacobian IS the derivative), so we only need
# a single FunctionWrapper variant.
function wrapfun_jac_iip(jac_f, inputs)
    return FunctionWrappersWrappers.FunctionWrappersWrapper(
        Void(jac_f), (typeof(inputs),), (Nothing,)
    )
end
