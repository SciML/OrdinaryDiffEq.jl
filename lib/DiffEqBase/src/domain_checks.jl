"""
    domain_checks_failing(checks, u, p, t)

Evaluate `checks` against `(u, p, t)` and return one message per failing check, naming each by
its position and by its source expression when it has one, or `nothing` if all pass. Every
check is evaluated rather than stopping at the first failure; the all-passing case allocates
nothing. Each check is `(u, p, t) -> Bool` or
`(u, p, t) -> Pair{Bool, <:AbstractString}`; a plain `Bool` uses a generic default failure
message. Used both by [`DomainCheckedFunction`](@ref) (to decide whether to call the real
right-hand-side) and by diagnostics (to explain which predicate failed).
"""
@inline function domain_checks_failing(checks, u, p, t)
    failures = nothing
    for (idx, check) in enumerate(checks)
        result = check(u, p, t)
        ok, msg = result isa Pair ? (result.first, result.second) : (result, "Domain check failed")
        ok && continue
        line = check isa TracedPredicate ? "[$idx] $msg: $(check.src)" : "[$idx] $msg"
        failures === nothing ? (failures = [line]) : push!(failures, line)
    end
    return failures
end

"""
    DomainCheckedFunction{iip, F, C}(f, pre_checks)

Wraps a raw ODE right-hand-side callable `f` together with a tuple/vector of `pre_checks`
predicates, each of the form `(u, p, t) -> Bool` or `(u, p, t) -> Pair{Bool, <:AbstractString}`.
Before calling `f`, every predicate is evaluated against `u`. If any predicate fails, `f` is
not called at all, and instead, a `NaN` is injected into the state. A failing check therefore 
falls into the NaN-detection in the adaptive step controller
(`integrator.isout = ... || isnan(EEst) || isinf(EEst)`) straight into a
shrink-and-retry.

Constructed automatically from the `domain_checks` solver keyword.
"""
struct DomainCheckedFunction{iip, F, C}
    f::F
    pre_checks::C
end

function DomainCheckedFunction{iip}(f::F, pre_checks::C) where {iip, F, C}
    return DomainCheckedFunction{iip, F, C}(f, pre_checks)
end

function (dcf::DomainCheckedFunction{true})(du, u, p, t)
    if domain_checks_failing(dcf.pre_checks, u, p, t) === nothing
        dcf.f(du, u, p, t)
    else
        du .= NaN
    end
    return nothing
end

function (dcf::DomainCheckedFunction{false})(u, p, t)
    domain_checks_failing(dcf.pre_checks, u, p, t) === nothing && return dcf.f(u, p, t)
    #maintain units of u
    return u .* (zero(eltype(u)) / zero(eltype(u)))
end


#part 2: deal with isoutofdomain for logging purposes
"""
    TracedPredicateLeaf(f, src, operand_srcs, operand_f)

One atomic sub-expression of an `isoutofdomain` predicate built by [`@isoutofdomain`](@ref),
retaining its source text so diagnostics can name it.

# Fields
- `f`: `(u, p, t) -> Bool` evaluating just this sub-expression.
- `src`: the sub-expression's source text, e.g. `"u[1] < 0"`.
- `operand_srcs`: source text of each operand, e.g. `("u[1]", "0")`. Empty when the leaf is
  not a comparison.
- `operand_f`: `(u, p, t) -> Tuple` of the operand *values*, or `nothing` when the leaf is
  not a comparison and so has no operands to report.
"""
struct TracedPredicateLeaf{F, S, O}
    f::F
    src::String
    operand_srcs::S
    operand_f::O
end

"""
    TracedPredicate(f, src, leaves)

An `isoutofdomain` predicate that carries its own source text, built by
[`@isoutofdomain`](@ref).

Purely a logging device: `f` is the user's original closure, stored verbatim and called
unchanged, so this costs nothing on the integrator's hot path (a single-field forward through
an immutable struct is erased by the compiler — the emitted code is instruction-for-instruction
identical to calling `f` directly). `src` and `leaves` are touched only by
[`isoutofdomain_report`](@ref) when a solve fails, so the tree never participates in deciding
whether a step is rejected.

It has to be the object stored in `opts.isoutofdomain` because that is the only handle
diagnostics have to reach the metadata from an integrator.
"""
struct TracedPredicate{F, L}
    f::F
    src::String
    leaves::L
end

@inline (tp::TracedPredicate)(u, p, t) = tp.f(u, p, t)

const COMPARISON_OPS = (:<, :>, :<=, :>=, :(==), :!=, :≤, :≥, :≠, :isless)

function is_connective(ex)
    ex isa Expr || return false
    (ex.head === :|| || ex.head === :&&) && return true
    ex.head === :call || return false
    length(ex.args) == 3 && ex.args[1] in (:|, :&) && return true
    return length(ex.args) == 2 && ex.args[1] === :!
end

# Strip LineNumberNodes and single-statement `begin ... end` wrappers, which the parser
# inserts around lambda bodies.
function strip_lines(ex)
    ex isa Expr || return ex
    if ex.head === :block
        stmts = filter(a -> !(a isa LineNumberNode), ex.args)
        length(stmts) == 1 && return strip_lines(only(stmts))
    end
    return Expr(ex.head, map(strip_lines, ex.args)...)
end

# Flatten a boolean expression into its atomic sub-expressions, in source order.
function collect_leaves!(acc, ex)
    if is_connective(ex)
        operands = (ex.head === :call) ? ex.args[2:end] : ex.args
        for operand in operands
            collect_leaves!(acc, operand)
        end
    else
        push!(acc, ex)
    end
    return acc
end

# Operands of a comparison leaf worth reporting the value of, or `nothing` if the leaf isn't a
# comparison or has no such operands.
function leaf_operands(ex)
    ex isa Expr || return nothing
    operands = if ex.head === :call && length(ex.args) == 3 && ex.args[1] in COMPARISON_OPS
        (ex.args[2], ex.args[3])
    elseif ex.head === :comparison
        # `0 < u[1] < 1` => (0, :<, u[1], :<, 1); operands sit at the odd indices
        Tuple(ex.args[1:2:end])
    else
        return nothing
    end
    reportable = filter(o -> !(o isa Union{Number, Bool, AbstractString, Char}), operands)
    return isempty(reportable) ? nothing : Tuple(reportable)
end

"""
    @isoutofdomain (u, p, t) -> expr

Build a [`TracedPredicate`](@ref) from an `isoutofdomain` predicate, retaining the source
text of the whole expression and of each atomic sub-expression so that a failed solve can
report *which* part of the predicate rejected the step, and with what operand values:

```julia
sol = solve(prob, Tsit5();
    isoutofdomain = @isoutofdomain (u, p, t) -> u[1] < 0 || u[2] > 10)
```

On a `dt <= dtmin` abort the diagnostics then read:

```
predicate: u[1] < 0 || u[2] > 10
  [1] u[1] < 0  => false  (u[1] = 0.42)
  [2] u[2] > 10 => true   (u[2] = 14.7)
```

The predicate is stored and called verbatim, so this is free at run time; the retained source
is read only when a solve fails. 

The argument must be an anonymous function of three parameters; their names are taken from
the signature, so `(state, par, time) -> state[1] < 0` reports `state[1] < 0`. Collection
recurses through `||`, `&&`, `!`, `|` and `&`; anything else (e.g. `any(x -> x < 0, u)`)
becomes a single leaf, still labelled with its source text but without operand decomposition.
"""
macro isoutofdomain(ex)
    if !(ex isa Expr && ex.head === :->)
        throw(ArgumentError("@isoutofdomain expects an anonymous function, e.g. `@isoutofdomain (u, p, t) -> u[1] < 0`; got `$ex`"))
    end
    signature = ex.args[1]
    argnames = signature isa Expr && signature.head === :tuple ? signature.args : [signature]
    if length(argnames) != 3 || !all(a -> a isa Symbol, argnames)
        throw(ArgumentError("@isoutofdomain expects exactly three plain parameters `(u, p, t)`; got `$signature`"))
    end

    body = strip_lines(ex.args[2])
    argtuple = Expr(:tuple, argnames...)

    leaf_exprs = map(collect_leaves!(Any[], body)) do leaf
        operands = leaf_operands(leaf)
        operand_srcs = operands === nothing ? () : Tuple(string.(operands))
        operand_f = operands === nothing ? nothing :
            esc(Expr(:->, argtuple, Expr(:tuple, operands...)))
        # `esc` so the sub-expression's free variables resolve in the caller's scope: a
        # predicate like `u[1] < lo` must find the caller's `lo`, not DiffEqBase's.
        leaf_f = esc(Expr(:->, argtuple, leaf))
        return :($TracedPredicateLeaf($leaf_f, $(string(leaf)), $operand_srcs, $operand_f))
    end

    return :($TracedPredicate($(esc(ex)), $(string(body)), ($(leaf_exprs...),)))
end

format_operand(v::Real) = @sprintf("%.4g", v)
format_operand(v) = repr(v)

function explain_leaf(idx, leaf, u, p, t, width)
    label = "  [$idx] " * rpad(leaf.src, width) * " => "
    value = try
        leaf.f(u, p, t)
    catch err
        return label * "errored: " * string(typeof(err))
    end
    operands = ""
    if leaf.operand_f !== nothing
        operands = try
            values = leaf.operand_f(u, p, t)
            pairs = ("$src = $(format_operand(val))" for (src, val) in zip(leaf.operand_srcs, values))
            "  (" * join(pairs, ", ") * ")"
        catch err
            "  (operands unavailable: " * string(typeof(err)) * ")"
        end
    end
    return label * string(value) * operands
end

"""
    isoutofdomain_report(pred, u, p, t) -> Vector{String}

Explain what an `isoutofdomain` predicate does at the state `(u, p, t)`,
for a diagnostic message. Returns the predicate's source text followed by one line per
atomic sub-expression giving its value and operand values.
"""
function isoutofdomain_report(tp::TracedPredicate, u, p, t)
    lines = ["predicate: " * tp.src]
    isempty(tp.leaves) && return lines
    width = maximum(length(leaf.src) for leaf in tp.leaves)
    for (idx, leaf) in enumerate(tp.leaves)
        push!(lines, explain_leaf(idx, leaf, u, p, t, width))
    end
    return lines
end

isoutofdomain_report(::Any, ::Any, ::Any, ::Any) = String[]
