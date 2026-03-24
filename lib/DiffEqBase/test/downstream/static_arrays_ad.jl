using OrdinaryDiffEq, StaticArrays, Test
using SciMLSensitivity  # Required for reverse-mode AD
using DifferentiationInterface
using ADTypes: AutoForwardDiff, AutoMooncake

# Load backends for all versions (required for DifferentiationInterface extensions)
using ForwardDiff
using Mooncake

# Version-dependent imports
if VERSION <= v"1.11"
    using Zygote
    using ADTypes: AutoZygote
end
if VERSION <= v"1.11"
    using Enzyme
    using ADTypes: AutoEnzyme
end

# Define backends based on Julia version
# ForwardDiff: All versions
# Mooncake: All versions
# Zygote: Julia <= 1.11
# Enzyme: Julia <= 1.11
function get_test_backends()
    backends = Pair{String, Any}[]
    # ForwardDiff on all versions
    push!(backends, "ForwardDiff" => AutoForwardDiff())
    # Mooncake on all versions
    push!(backends, "Mooncake" => AutoMooncake(; config = nothing))
    # Zygote only on Julia <= 1.11
    if VERSION <= v"1.11"
        push!(backends, "Zygote" => AutoZygote())
    end
    # Enzyme only on Julia <= 1.11
    if VERSION <= v"1.11"
        push!(backends, "Enzyme" => AutoEnzyme())
    end
    return backends
end

f(u, p, t) = copy(u)

backends = get_test_backends()
# Note: Mooncake and Zygote are excluded due to issues with SciMLSensitivity adjoint rules
# Mooncake: issues with adjoint rules
# Zygote: MethodError length(::Nothing) in automatic_sensealg_choice when p=nothing
backends_for_staticarrays = filter(b -> b[1] ∉ ("Mooncake", "Zygote"), backends)

@testset "StaticArrays AD tests" begin
    for (name, backend) in backends_for_staticarrays
        @testset "StaticArrays derivative with $name" begin
            du1 = DifferentiationInterface.derivative(backend, 5.0) do x
                prob = ODEProblem(f, [x], (0.0, 1.0), nothing)
                sol = solve(prob, Tsit5())
                sol.u[end][1]
            end

            du2 = DifferentiationInterface.derivative(backend, 5.0) do x
                prob = ODEProblem(f, SVector(x), (0.0, 1.0), nothing)
                sol = solve(prob, Tsit5())
                sol.u[end][1]
            end

            @test du1 ≈ du2
        end
    end
end
