struct NoIndexArray{T, N} <: AbstractArray{T, N}
    x::Array{T, N}
end
Base.size(x::NoIndexArray) = size(x.x)
Base.axes(x::NoIndexArray) = axes(x.x)
function Base.similar(x::NoIndexArray, dims::Union{Integer, AbstractUnitRange}...)
    NoIndexArray(similar(x.x, dims...))
end
Base.copyto!(x::NoIndexArray, y::NoIndexArray) = NoIndexArray(copyto!(x.x, y.x))
Base.copy(x::NoIndexArray) = NoIndexArray(copy(x.x))
Base.zero(x::NoIndexArray) = NoIndexArray(zero(x.x))
Base.fill!(x::NoIndexArray, y) = NoIndexArray(fill!(x.x, y))
Base.mapreduce(f, op, x::NoIndexArray; kwargs...) = mapreduce(f, op, x.x; kwargs...)
Base.any(f::Function, x::NoIndexArray; kwargs...) = any(f, x.x; kwargs...)
Base.all(f::Function, x::NoIndexArray; kwargs...) = all(f, x.x; kwargs...)
Base.:(==)(x::NoIndexArray, y::NoIndexArray) = x.x == y.x

struct NoIndexStyle{N} <: Broadcast.AbstractArrayStyle{N} end
NoIndexStyle(::Val{N}) where {N} = NoIndexStyle{N}()
NoIndexStyle{M}(::Val{N}) where {N, M} = NoIndexStyle{N}()
Base.BroadcastStyle(::Type{<:NoIndexArray{T, N}}) where {T, N} = NoIndexStyle{N}()
function Base.similar(bc::Base.Broadcast.Broadcasted{NoIndexStyle{N}},
    ::Type{ElType}) where {N, ElType}
    NoIndexArray(similar(Array{ElType, N}, axes(bc)))
end
Base.Broadcast._broadcast_getindex(x::NoIndexArray, i) = x.x[i]
Base.Broadcast.extrude(x::NoIndexArray) = x
using ArrayInterface
ArrayInterface.fast_scalar_indexing(::Type{<:NoIndexArray}) = false

@inline function Base.copyto!(dest::NoIndexArray,
    bc::Base.Broadcast.Broadcasted{<:NoIndexStyle})
    axes(dest) == axes(bc) || throwdm(axes(dest), axes(bc))
    bc′ = Base.Broadcast.preprocess(dest, bc)
    dest′ = dest.x
    @simd for I in eachindex(bc′)
        @inbounds dest′[I] = bc′[I]
    end
    return dest
end
Base.show_vector(io::IO, x::NoIndexArray) = Base.show_vector(io, x.x)

Base.show(io::IO, x::NoIndexArray) = (print(io, "NoIndexArray"); show(io, x.x))
function Base.show(io::IO, ::MIME"text/plain", x::NoIndexArray)
    println(io, Base.summary(x), ":")
    Base.print_array(io, x.x)
end

using OrdinaryDiffEq
prob = ODEProblem((du, u, p, t) -> copyto!(du, u), NoIndexArray(ones(10, 10)), (0.0, 10.0))
algs = [Tsit5(), BS3(), Vern9(), DP5()]
for alg in algs
    sol = @test_nowarn solve(prob, alg)
    @test_nowarn sol(0.1)
    @test_nowarn sol(similar(prob.u0), 0.1)
end
