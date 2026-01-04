alg_order(::ExplicitTaylor2) = 2
alg_stability_size(alg::ExplicitTaylor2) = 1

alg_order(::ExplicitTaylor{P}) where {P} = P
alg_stability_size(alg::ExplicitTaylor) = 1

JET_CACHE = IdDict()

function build_jet(f::ODEFunction{iip}, p, order, length = nothing) where {iip}
    return build_jet(f, Val{iip}(), p, order, length)
end

function build_jet(f, ::Val{iip}, p, order::Val{P}, length = nothing) where {P, iip}
    if haskey(JET_CACHE, f)
        list = JET_CACHE[f]
        index = findfirst(x -> x[1] == order && x[2] == p, list)
        index !== nothing && return list[index][3]
    end
    @variables t0::Real
    u0 = isnothing(length) ? Symbolics.variable(:u0) : Symbolics.variables(:u0, 1:length)
    if iip
        @assert length isa Integer
        f0 = similar(u0)
        f(f0, u0, p, t0)
    else
        f0 = f(u0, p, t0)
    end
    u = TaylorDiff.make_seed(u0, f0, Val(1))
    for index in 2:P
        t = TaylorScalar{index - 1}(t0, one(t0))
        if iip
            fu = similar(u)
            f(fu, u, p, t)
        else
            fu = f(u, p, t)
        end
        d = get_coefficient(fu, index - 1) / index
        u = append_coefficient(u, d)
    end
    jet = build_function(u, u0, t0; expression = Val(false), cse = true)
    if !haskey(JET_CACHE, f)
        JET_CACHE[f] = []
    end
    push!(JET_CACHE[f], (order, p, jet))
    return jet
end

# evaluate using Qin Jiushao's algorithm
@generated function evaluate_polynomial(t::TaylorScalar{T, P}, z) where {T, P}
    ex = :(v[$(P + 1)])
    for i in P:-1:1
        ex = :(v[$i] + z * $ex)
    end
    return :($(Expr(:meta, :inline)); v = flatten(t); $ex)
end
