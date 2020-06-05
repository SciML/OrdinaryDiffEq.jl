@kernel function arrayfuse_linear(af, src)
	i = @index(Global)
	af.visible[i] = af.p[1] * af.visible[i] + af.p[2] * src[i]
	af.hidden[i] = af.hidden[i] + af.p[3] * af.visible[i]
end

get_device(::Type) = CPU()
# TODO: Add dispatch for GPU
# get_device(::Type{<:CuArray}) = GPU()
# but should we add CuArray to dependency?

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
struct ArrayFuse{AT, T, P} <: AbstractArray{T, 1}
	visible::AT
	hidden::AT
	p::P
end

ArrayFuse(visible::AT, hidden::AT, p) where {AT} = ArrayFuse{AT, eltype(visible), typeof(p)}(visible, hidden, p)

@inline function Base.copyto!(af::OrdinaryDiffEq.ArrayFuse{AT, T, P}, src::Base.Broadcast.Broadcasted) where {AT, T, P}
	device = get_device(AT)
	kernel = arrayfuse_linear(device, 16)
	event = kernel(af, src; ndrange=length(af.visible))
	wait(device, event)
end

# not recommended but good to have
@inline function Base.getindex(af::ArrayFuse, index)
	return af.visible[index]
end

@inline function Base.setindex!(af::ArrayFuse, value, index)
	af.visible[index] = af.p[1] * af.visible[index] + af.p[2] * value
	af.hidden[index] =  muladd(af.p[3], af.visible[index], af.hidden[index])
end

@inline Base.size(af::ArrayFuse) = length(af.visible)
@inline Base.axes(af::ArrayFuse) = axes(af.visible)

Adapt.adapt_structure(to, af::ArrayFuse{AT, T, P}) where {AT, T, P} = 
	ArrayFuse(adapt(to, af.visible), adapt(to, af.hidden), af.p)
