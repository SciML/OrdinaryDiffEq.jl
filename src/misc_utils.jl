immutable DiffCache{T, S}
    du::Vector{T}
    dual_du::Vector{S}
end

function DiffCache{chunk_size}(T, length, ::Type{Val{chunk_size}})
    DiffCache(zeros(T, length), zeros(Dual{chunk_size, T}, length))
end

DiffCache(u::AbstractArray) = DiffCache(eltype(u),length(u),Val{ForwardDiff.pickchunksize(u)})

get_du{T<:Dual}(dc::DiffCache, ::Type{T}) = dc.dual_du
get_du(dc::DiffCache, T) = dc.du

type MutableReference{T}
val::T
end

#Base.getindex(m::MutableReference)=m.val
#Base.setindex!(m::MutableReference{T}, v::{T},I...) = m.val=v
