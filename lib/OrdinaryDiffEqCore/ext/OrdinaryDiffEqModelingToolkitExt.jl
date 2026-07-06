module OrdinaryDiffEqModelingToolkitExt

using OrdinaryDiffEqCore, ModelingToolkit
using Printf: @sprintf

function OrdinaryDiffEqCore.system_singularity_rootcause(sys, u)
    substitution_map = Dict(zip(unknowns(sys), u))
    diagnosis = String[]
    for eq in equations(sys)
        find_singular_subterms(eq, eq.rhs, substitution_map, diagnosis)
    end
    return diagnosis
end

function find_singular_subterms(eq, expr, sub_map, diagnosis)
    expr = Symbolics.unwrap(expr)
    !SymbolicUtils.iscall(expr) && return diagnosis
    op = SymbolicUtils.operation(expr)
    args = SymbolicUtils.arguments(expr)

    if op === (/) #division, singular if we divide by small thing
        d = Symbolics.value(Symbolics.substitute(args[2], sub_map))
        if d isa Number && abs(d) < 1e-10
            push!(diagnosis, "in equation $eq: division by very small value $(args[2]) ≈ $(@sprintf("%.4g", d)) leads to singularity.")
        end
    elseif op === log #singular if we log small thing
        x = Symbolics.value(Symbolics.substitute(args[1], sub_map))
        if x isa Number && x <= 1e-10
            push!(diagnosis, "in equation $eq: log of $(args[1]) = $(@sprintf("%.4g", x)) near/at singularity (derivative blows up).")
        end
    elseif op === sqrt 
        x = Symbolics.value(Symbolics.substitute(args[1], sub_map))
        if x isa Number && x < 1e-10
            push!(diagnosis, "in equation $eq: sqrt of $(args[1]) = $(@sprintf("%.4g", x)) near/at singularity (derivative blows up).")
        end
    elseif op === (^)
        e = Symbolics.value(Symbolics.substitute(args[2], sub_map))
        b = Symbolics.value(Symbolics.substitute(args[1], sub_map))
        if e isa Number && b isa Number #two cases
            if e < 0 && abs(b) < 1e-10
                push!(diagnosis, "in equation $eq: ($(args[1])) raised to power $e with base ≈ $(@sprintf("%.4g", b)) going to 0; result diverges.")
            elseif e > 0 && abs(b) > 1
                push!(diagnosis, "in equation $eq: ($(args[1]) ≈ $(@sprintf("%.4g", b))) raised to power $e - base magnitude is large and being amplified.")
            end
        end
    end

    for arg in args
        find_singular_subterms(eq, arg, sub_map, diagnosis)
    end
    return diagnosis
end

end
