"""
    ArrayFuse{AT, T, P} <: AbstractArray{T, 1}

GPU Friendly type to wrap around two arrays - `visible` and `hidden`, for which when we `setindex!` some value `v` at index `i`
we get

```
visible[i] = p[1] * visible[i] + p[2] * v
hidden[i] = hidden[i] + p[3] * visible[i]
```

where p is a parameter tuple of size 3.
"""
struct ArrayFuse{AT, T, P}
    visible::AT
    hidden::AT
    p::P
end

Base.ndims(::Type{<:ArrayFuse{AT}}) where {AT} = ndims(AT)

function ArrayFuse(visible::AT, hidden::AT, p) where {AT}
    ArrayFuse{AT, eltype(visible), typeof(p)}(visible, hidden, p)
end

function Base.copyto!(af::ArrayFuse, src::Broadcast.Broadcasted)
    @. af.visible = af.p[1] * af.visible + af.p[2] * src
    @. af.hidden = af.hidden + af.p[3] * af.visible
end

function Base.materialize!(af::ArrayFuse, src::Broadcast.Broadcasted)
    copyto!(af, src)
end

# not recommended but good to have
function Base.getindex(af::ArrayFuse, index::Int)
    return af.visible[index]
end

function Base.setindex!(af::ArrayFuse, value, index::Int)
    af.visible[index] = af.p[1] * af.visible[index] + af.p[2] * value
    af.hidden[index] = muladd(af.p[3], af.visible[index], af.hidden[index])
end

Base.size(af::ArrayFuse) = length(af.visible)
Base.axes(af::ArrayFuse) = axes(af.visible)

function Adapt.adapt_structure(to, af::ArrayFuse{AT, T, P}) where {AT, T, P}
    ArrayFuse(adapt(to, af.visible), adapt(to, af.hidden), af.p)
end
