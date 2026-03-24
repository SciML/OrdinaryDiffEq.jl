using LinearAlgebra, OrdinaryDiffEq, Test
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

# setup
pd = 3

## ode with complex numbers
H0 = rand(ComplexF64, pd, pd)
A = rand(ComplexF64, pd, pd)
function f!(du, u, p, t)
    a, b, c = p
    du .= (A * u) * (a * cos(b * t + c))
    du .+= H0 * u
    return nothing
end

## time span
tspan = (0.0, 1.0)

## initial state
u0 = hcat(normalize(rand(ComplexF64, pd)), normalize(rand(pd)))

## ode problem
prob0 = ODEProblem(
    f!, u0, tspan, rand(3); saveat = range(tspan..., length = 3),
    reltol = 1.0e-6,
    alg = Tsit5()
)
## final state cost
cost(u) = abs2(tr(first(u)'u[2])) - abs2(tr(first(u)'last(u)))

## real loss function via complex ode
function loss(p)
    prob = remake(prob0; p)
    sol = solve(prob)
    return cost(sol.u) + sum(p) / 10
end

## same problem via reals
### realify complex ode problem
function real_f(du, u, p, t)
    complex_u = complex.(selectdim(u, 3, 1), selectdim(u, 3, 2))
    complex_du = copy(complex_u)
    prob0.f(complex_du, complex_u, p, t)
    selectdim(du, 3, 1) .= real(complex_du)
    selectdim(du, 3, 2) .= imag(complex_du)
    return nothing
end
prob0_real = remake(prob0; f = real_f, u0 = cat(real(prob0.u0), imag(prob0.u0); dims = 3))
### real loss function via real ode
function loss_via_real(p)
    prob = remake(prob0_real; p)
    sol = solve(prob)
    u = [complex.(selectdim(u, 3, 1), selectdim(u, 3, 2)) for u in sol.u]
    return cost(u) + sum(p) / 10
end

# assert
@assert eltype(last(solve(prob0).u)) <: Complex
@assert eltype(last(solve(prob0_real).u)) <: Real
function assert_fun()
    p0 = rand(3)
    return isapprox(loss(p0), loss_via_real(p0); rtol = 1.0e-4)
end
@assert all([assert_fun() for _ in 1:(2^6)])

# test ad with DifferentiationInterface using multiple backends
backends = get_test_backends()
# Note: Mooncake and Enzyme are excluded due to issues with complex number ODE differentiation
backends_for_complex = filter(b -> b[1] ∉ ("Mooncake", "Enzyme"), backends)

function test_ad_with_backend(backend, name)
    p0 = rand(3)
    grad_real = DifferentiationInterface.gradient(loss_via_real, backend, p0)
    grad_complex = DifferentiationInterface.gradient(loss, backend, p0)
    any(isnan.(grad_complex)) &&
        @warn "NaN detected in gradient using ode with complex numbers with $name !!"
    any(isnan.(grad_real)) && @warn "NaN detected in gradient using realified ode with $name !!"
    rel_err = norm(grad_complex - grad_real) / max(norm(grad_complex), norm(grad_real))
    return isapprox(grad_complex, grad_real; rtol = 1.0e-6) ? true : (@show name, rel_err; false)
end

@testset "Complex number AD tests" begin
    for (name, backend) in backends_for_complex
        @testset "AD via ode with complex numbers ($name)" begin
            @time @test all([test_ad_with_backend(backend, name) for _ in 1:(2^6)])
        end
    end
end
