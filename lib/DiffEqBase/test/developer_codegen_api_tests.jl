using DiffEqBase: DiffEqBase, @tight_loop_macros
using Test

module ExternalDiffEqBaseCodegen
    using DiffEqBase: @tight_loop_macros

    function add_one!(out, x)
        @tight_loop_macros for i in eachindex(out, x)
            @inbounds out[i] = x[i] + 1
        end
        return out
    end
end

@test ExternalDiffEqBaseCodegen.add_one!(zeros(3), [1.0, 2.0, 3.0]) == [2.0, 3.0, 4.0]

@static if isdefined(Base, :ispublic)
    @test Base.ispublic(DiffEqBase, Symbol("@tight_loop_macros"))
    @test Base.Docs.hasdoc(DiffEqBase, Symbol("@tight_loop_macros"))
end
