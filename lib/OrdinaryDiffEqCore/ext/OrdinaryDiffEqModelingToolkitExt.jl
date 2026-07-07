module OrdinaryDiffEqModelingToolkitExt

using OrdinaryDiffEqCore, ModelingToolkit
using Printf: @sprintf

function OrdinaryDiffEqCore.system_singularity_rootcause(sys, u, uprev)
    prev_substitution_map = Dict(zip(unknowns(sys), uprev))
    diagnosis = String[]
    for eq in equations(sys)
        find_singular_subterms(eq, eq.rhs, prev_substitution_map, diagnosis)
    end

    #check for assertion failures
    unks = unknowns(sys)
    curr_substitution_map = Dict(zip(unks, u))
    for (cond, msg) in ModelingToolkit.assertions(sys)
        f = Symbolics.build_function(cond, unks; expression = Val{false})
        if f(u) === false
            push!(diagnosis, "\nAssertion violated: $cond - \"$msg\"")
            find_failing_subterms(cond, prev_substitution_map, curr_substitution_map, diagnosis)
        end
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

function find_failing_subterms(cond, prev_map, curr_map, diagnosis)
    c = Symbolics.unwrap(cond)
    !SymbolicUtils.iscall(c) && return diagnosis
    op = SymbolicUtils.operation(c)
    args = SymbolicUtils.arguments(c)

    if (op === (<) || op === (>) || op === (<=) || op === (>=)) && length(args) == 2
        #compare using previous non-nan values to find violating subclauses, then output current values
        lhs = Symbolics.value(Symbolics.substitute(args[1], prev_map))
        rhs = Symbolics.value(Symbolics.substitute(args[2], prev_map))
        if lhs isa Number && rhs isa Number
            # small margin -> violated
            margin = (op === (<) || op === (<=)) ? rhs - lhs : lhs - rhs
            if margin <= 1e-6
                push!(diagnosis, "   subclause `$c` violated: $(clause_values(c, curr_map))")
            end
        end
    elseif op === (!=) && length(args) == 2
        lhs = Symbolics.value(Symbolics.substitute(args[1], prev_map))
        rhs = Symbolics.value(Symbolics.substitute(args[2], prev_map))
        if lhs isa Number && rhs isa Number && abs(lhs - rhs) <= 1e-6
            push!(diagnosis, "   subclause `$c` violated: $(clause_values(c, curr_map))")
        end
    else #recurse
        for arg in args
            find_failing_subterms(arg, prev_map, curr_map, diagnosis)
        end
    end
    return diagnosis
end

function clause_values(c, curr_map)
    parts = String[]
    for v in Symbolics.get_variables(c)
        val = Symbolics.value(Symbolics.substitute(v, curr_map))
        push!(parts, val isa Number ? "$v = $(@sprintf("%.4g", val))" : "$v = $val")
    end
    return join(parts, ", ")
end

end
