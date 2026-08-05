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

Constructed automatically from the `domain_checks` solver keyword, by
[`apply_domain_checks`](@ref), as the outermost layer of `f.f` — outside any
`FunctionWrappersWrapper` `promote_f` may have installed. [`strip_domain_checks`](@ref)
removes it again for the code paths that must not see the checks.
"""
struct DomainCheckedFunction{iip, F, C}
    f::F
    pre_checks::C
end

function DomainCheckedFunction{iip}(f::F, pre_checks::C) where {iip, F, C}
    return DomainCheckedFunction{iip, F, C}(f, pre_checks)
end

unwrapped_f(dcf::DomainCheckedFunction) = unwrapped_f(dcf.f)

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

"""
    strip_domain_checks(f)

Remove a [`DomainCheckedFunction`](@ref) layer, returning `f` unchanged when there is none.
Accepts either the inner callable or the `AbstractSciMLFunction` wrapping it.

Checks gate *states* the integrator might accept — RK stages, Newton residuals — which all
go through `integrator.f` and keep them. They must not gate *derivative probes*: the point
`uprev ± ε eᵢ` a finite-difference Jacobian evaluates is not a state, `uprev` itself is
in-domain, and only the stencil crossed the boundary. Poisoning that Jacobian column would
reject a valid step unrecoverably, since `ε ≈ sqrt(eps) * |u|` does not scale with `dt`, so
every retry probes the same point and the solve marches to `dtmin`, and only under
`AutoFiniteDiff`, since `AutoForwardDiff` perturbs the dual part and leaves every predicate
comparison in-domain.

So the differentiation wrappers (`UJacobianWrapper`/`UDerivativeWrapper`/`TimeGradientWrapper`
and their `jac_config`s) are built from a stripped `f`, and `islinearfunction` strips before
asking `islinear`, since it hands `f.f` to the solver *as* the Jacobian.
"""
strip_domain_checks(dcf::DomainCheckedFunction) = dcf.f
strip_domain_checks(f) = f
function strip_domain_checks(f::SciMLBase.AbstractSciMLFunction)
    hasfield(typeof(f), :f) || return f
    inner = strip_domain_checks(f.f)
    inner === f.f && return f
    return @set f.f = inner
end

"""
    find_domain_checks(f) -> checks or nothing

Recover the `domain_checks` predicates carried by `f`, or `nothing` if it carries none.
Used by failure diagnostics, which have only the integrator to work from. Kept separate
from `f.f isa DomainCheckedFunction` so that the diagnostics do not silently go quiet if
the layering around `f.f` ever changes.
"""
find_domain_checks(dcf::DomainCheckedFunction) = dcf.pre_checks
find_domain_checks(::Any) = nothing
function find_domain_checks(f::SciMLBase.AbstractSciMLFunction)
    # See `strip_domain_checks`: not every `AbstractSciMLFunction` has a single `f`.
    hasfield(typeof(f), :f) || return nothing
    return find_domain_checks(f.f)
end

"""
    supports_domain_checks(alg) -> Bool

Whether `alg` implements the step-rejection recovery that the `domain_checks` keyword
relies on. `false` by default: DiffEqBase installs the checks during problem
concretization, which every solver package routes through, but only the solvers that turn
a NaN right-hand-side into a smaller-`dt` retry can honor them. Passing `domain_checks` to
any other solver is an error rather than a silently different meaning.
"""
supports_domain_checks(alg) = false

"""
    apply_domain_checks(f, domain_checks, alg; opaque_params = false) -> f

Install the `domain_checks` predicates on `f` as the outermost layer of `f.f`, or return
`f` unchanged when no checks were requested. Called from `get_concrete_problem`, i.e.
before the cache, the nonlinear solver and the solution object are built, so that every
consumer of `f` sees one consistent function.

Rejects the combinations that cannot honor the checks instead of silently dropping them:
right-hand-sides whose call signature is not `(du, u, p, t)`/`(u, p, t)` (split, DAE and
second-order/partitioned forms, which would `MethodError`, or evaluate only one of several
sub-functions), operator right-hand-sides (where the solver uses `f.f` itself as the
Jacobian), and algorithms without the retry path.
"""
function apply_domain_checks(f, domain_checks, alg; opaque_params::Bool = false)
    (domain_checks === nothing || isempty(domain_checks)) && return f

    if !supports_domain_checks(alg)
        throw(
            ArgumentError(
                "domain_checks is not supported for $(alg === nothing ? "this solver" : nameof(typeof(alg))). " *
                    "The keyword relies on a failing predicate being converted into a step retry with a " *
                    "smaller dt, which only the OrdinaryDiffEq.jl solvers implement."
            )
        )
    end
    if !(f isa ODEFunction)
        throw(
            ArgumentError(
                "domain_checks is only supported for ODEFunction right-hand-sides, got " *
                    "$(nameof(typeof(f))). Split (IMEX), DAE, and second-order/partitioned " *
                    "right-hand-sides either take a different call signature or dispatch to several " *
                    "sub-functions, so a single wrapped callable would miss evaluations rather than " *
                    "gate them."
            )
        )
    end

    if f.f isa SciMLBase.AbstractSciMLOperator || islinear(f.f)
        throw(
            ArgumentError(
                "domain_checks is not supported for operator or linear right-hand-sides. For a " *
                    "linear `f.f` the solvers use it directly as the Jacobian (see " *
                    "`islinearfunction`), which wrapping would silently replace with a " *
                    "numerically-formed dense Jacobian."
            )
        )
    end
    if opaque_params
        throw(
            ArgumentError(
                "domain_checks is not supported together with AutoDePSpecialize's opaque parameter " *
                    "packing: the predicates sit outside the function wrapper that unpacks `p`, so they " *
                    "would receive the opaque container rather than the problem's parameters."
            )
        )
    end
    # Idempotent: `get_concrete_problem` can run more than once on the same problem (e.g. a
    # sensitivity adjoint re-concretizing), and the checks must not stack up.
    find_domain_checks(f) === nothing || return f

    return @set f.f = DomainCheckedFunction{isinplace(f)}(f.f, domain_checks)
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
