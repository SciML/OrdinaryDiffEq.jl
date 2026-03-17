"""
    Discontinuity(t, order)

Object of discontinuity of order `order` at time `t`, i.e. discontinuity of `order`th
derivative at time `t`.
"""
struct Discontinuity{tType, O}
    t::tType
    order::O
end

# ordering of discontinuities
Base.:<(a::Discontinuity, b::Discontinuity) = a.t < b.t || (a.t == b.t && a.order < b.order)
function Base.isless(a::Discontinuity, b::Discontinuity)
    return isless(a.t, b.t) || (isequal(a.t, b.t) && isless(a.order, b.order))
end
function Base.isequal(a::Discontinuity, b::Discontinuity)
    return isequal(a.t, b.t) && isequal(a.order, b.order)
end
Base.:(==)(a::Discontinuity, b::Discontinuity) = a.t == b.t && a.order == b.order

# ordering with numbers
Base.:<(a::Discontinuity, b::Number) = a.t < b
Base.:<(a::Number, b::Discontinuity) = a < b.t
Base.isless(a::Discontinuity, b::Number) = isless(a.t, b)
Base.isless(a::Number, b::Discontinuity) = isless(a, b.t)
Base.:(==)(a::Discontinuity, b::Number) = a.t == b
Base.:(==)(a::Number, b::Discontinuity) = a == b.t
Base.isequal(a::Discontinuity, b::Number) = isequal(a.t, b)
Base.isequal(a::Number, b::Discontinuity) = isequal(a, b.t)

# multiplication with numbers
Base.:*(a::Discontinuity, b::Number) = Discontinuity(a.t * b, a.order)
Base.:*(a::Number, b::Discontinuity) = Discontinuity(a * b.t, b.order)

# conversion to numbers
Base.convert(::Type{T}, d::Discontinuity) where {T <: Number} = convert(T, d.t)
Base.convert(::Type{T}, d::Discontinuity{T}) where {T <: Number} = d.t

# Fix ForwardDiff ambiguities
function Base.convert(
        ::Type{ForwardDiff.Dual{T, V, N}},
        ::Discontinuity{ForwardDiff.Dual{T, V, N}}
    ) where {N, V, T}
    throw(
        MethodError(
            convert, (
                ForwardDiff.Dual{T, V, N}, Discontinuity{ForwardDiff.Dual{T, V, N}},
            )
        )
    )
end
function Base.convert(::Type{ForwardDiff.Dual{T, V, N}}, ::Discontinuity) where {N, V, T}
    throw(MethodError(convert, (ForwardDiff.Dual{T, V, N}, Discontinuity)))
end

# conversion to discontinuity
function Base.convert(::Type{Discontinuity{T}}, x::Number) where {T <: Number}
    return convert(Discontinuity{T, Int}, x)
end
function Base.convert(::Type{Discontinuity{T, O}}, x::Number) where {T <: Number, O}
    return Discontinuity{T, O}(convert(T, x), zero(O))
end
function Base.convert(::Type{Discontinuity{T, O}}, x::T) where {T <: Number, O}
    return Discontinuity{T, O}(x, zero(O))
end

# display
Base.show(io::IO, d::Discontinuity) = print(io, d.t, " (order ", d.order, ")")
Base.show(io::IO, ::MIME"text/plain", d::Discontinuity) = print(io, typeof(d), ":\n   ", d)
