using SciMLTesting

activate_group_env(
    @__DIR__;
    parent = [joinpath(@__DIR__, "..", ".."), joinpath(@__DIR__, "..", "..", "..", "..")],
)
include(joinpath(@__DIR__, "qa.jl"))
