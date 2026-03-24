using Test
using DiffEqBase
using Aqua

@testset "Aqua tests (performance)" begin
    # This tests that we don't accidentally run into
    # https://github.com/JuliaLang/julia/issues/29393
    # Aqua.test_unbound_args(DiffEqBase) # fails
    ua = Aqua.detect_unbound_args_recursively(DiffEqBase)
    @test length(ua) == 0
    # Uncomment for debugging:
    # @show ua

    # See: https://github.com/SciML/OrdinaryDiffEq.jl/issues/1750
    # Test that we're not introducing method ambiguities across deps
    ambs = Aqua.detect_ambiguities(DiffEqBase; recursive = true)
    pkg_match(pkgname, pkdir::Nothing) = false
    pkg_match(pkgname, pkdir::AbstractString) = occursin(pkgname, pkdir)
    filter!(x -> pkg_match("DiffEqBase", pkgdir(last(x).module)), ambs)

    # Uncomment for debugging:
    # for method_ambiguity in ambs
    #     @show method_ambiguity
    # end
    @warn "Number of method ambiguities: $(length(ambs))"
    @test length(ambs) ≤ 4
end

@testset "Aqua tests (additional)" begin
    Aqua.test_undefined_exports(DiffEqBase)
    Aqua.test_stale_deps(DiffEqBase)
    Aqua.test_deps_compat(DiffEqBase)
    Aqua.test_project_extras(DiffEqBase)
    Aqua.test_project_toml_formatting(DiffEqBase)
    # Aqua.test_piracy(DiffEqBase) # failing
end

nothing
