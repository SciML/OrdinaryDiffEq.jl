using OrdinaryDiffEqCore: OrdinaryDiffEqCore, @fold, @OnDemandTableauExtract
using Test

module ExternalOrdinaryDiffEqCodegen
    using OrdinaryDiffEqCore: @fold, @OnDemandTableauExtract

    struct OneTypeTableau{T}
        first::T
        second::T
    end
    OneTypeTableau(::Type{T}) where {T} = OneTypeTableau(one(T), convert(T, 1 // 2))

    struct TwoTypeTableau{T, T2}
        weight::T
        node::T2
    end
    TwoTypeTableau(::Type{T}, ::Type{T2}) where {T, T2} =
        TwoTypeTableau(one(T), convert(T2, 1 // 4))

    function one_type_coefficients(T)
        @OnDemandTableauExtract OneTypeTableau T
        return first, second
    end

    function two_type_coefficients(T, T2)
        @OnDemandTableauExtract TwoTypeTableau T T2
        return weight, node
    end

    @fold function folded_pair(::Type{T}) where {T}
        return convert(T, 1 // 2), one(T)
    end
end

@test ExternalOrdinaryDiffEqCodegen.one_type_coefficients(Float64) == (1.0, 0.5)
@test ExternalOrdinaryDiffEqCodegen.two_type_coefficients(Float64, Float32) ==
    (1.0, 0.25f0)
@test ExternalOrdinaryDiffEqCodegen.folded_pair(Float64) == (0.5, 1.0)

@static if isdefined(Base, :ispublic)
    for name in (Symbol("@fold"), Symbol("@OnDemandTableauExtract"))
        @test Base.ispublic(OrdinaryDiffEqCore, name)
        @test Base.Docs.hasdoc(OrdinaryDiffEqCore, name)
    end
end
