using Test

@testset "Sublibrary group environment precedence" begin
    runners = Dict(
        joinpath(@__DIR__, "..", "..", "lib", "DiffEqBase", "test", "runtests.jl") => "Core",
        joinpath(@__DIR__, "..", "..", "lib", "DiffEqDevTools", "test", "runtests.jl") => "ALL",
        joinpath(@__DIR__, "..", "..", "lib", "DelayDiffEq", "test", "runtests.jl") => "ALL",
        joinpath(@__DIR__, "..", "..", "lib", "StochasticDiffEq", "test", "runtests.jl") => "ALL",
    )

    for (runner, default) in runners
        assignment = match(r"const TEST_GROUP = ([^\n]+)", read(runner, String))
        @test assignment !== nothing
        selector = Meta.parse(only(assignment.captures))

        withenv("ODEDIFFEQ_TEST_GROUP" => nothing, "GROUP" => "Requested") do
            @test Core.eval(@__MODULE__, selector) == "Requested"
        end
        withenv("ODEDIFFEQ_TEST_GROUP" => "Specific", "GROUP" => "Requested") do
            @test Core.eval(@__MODULE__, selector) == "Specific"
        end
        withenv("ODEDIFFEQ_TEST_GROUP" => nothing, "GROUP" => nothing) do
            @test Core.eval(@__MODULE__, selector) == default
        end
    end
end
