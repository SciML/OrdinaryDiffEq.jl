struct CNABTableau{T,T2}
    a0::T
    a1::T
end

Base.@pure function CNABTableau{T<:CompiledFloats,T2<:CompiledFloats}(::Type{T},::Type{T2})
    a0 = T(-0.5)
    a1 = T(1.5) 
    CNABTableau(a0,a1)
end

function CNABTableau(T,T2)
    a0 = T(-1//2)
    a1 = T(3//2)
    CNABTableau(a0,a1)
end