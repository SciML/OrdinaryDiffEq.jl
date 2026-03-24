# Precompilation workload for DiffEqBase
# This precompiles commonly used code paths to reduce TTFX

PrecompileTools.@setup_workload begin
    # Minimal setup for precompilation
    u0 = [1.0, 0.0, 0.0]
    u1 = [1.1, 0.1, 0.1]
    t = 0.0
    dt = 0.1
    α = 1.0e-6
    ρ = 1.0e-3

    PrecompileTools.@compile_workload begin
        # Precompile ODE_DEFAULT_NORM for Vector{Float64} (most common case)
        ODE_DEFAULT_NORM(u0, t)

        # Precompile NAN_CHECK for Vector{Float64}
        NAN_CHECK(u0)

        # Precompile INFINITE_OR_GIANT for Vector{Float64}
        INFINITE_OR_GIANT(u0)

        # Precompile ODE_DEFAULT_UNSTABLE_CHECK for Vector{Float64}
        # Using nothing for p since that's very common
        ODE_DEFAULT_UNSTABLE_CHECK(dt, u0, nothing, t)

        # Precompile ODE_DEFAULT_ISOUTOFDOMAIN
        ODE_DEFAULT_ISOUTOFDOMAIN(u0, nothing, t)

        # Precompile ODE_DEFAULT_PROG_MESSAGE
        ODE_DEFAULT_PROG_MESSAGE(dt, u0, nothing, t)

        # Precompile recursive_length for Vector{Float64}
        recursive_length(u0)

        # Precompile UNITLESS_ABS2 for Vector{Float64}
        UNITLESS_ABS2(u0)

        # Precompile calculate_residuals for Vector{Float64} (most common case)
        calculate_residuals(u0, u1, α, ρ, ODE_DEFAULT_NORM, t)

        # Precompile in-place version
        out = similar(u0)
        calculate_residuals!(out, u0, u1, α, ρ, ODE_DEFAULT_NORM, t)

        # Precompile with explicit differences
        ũ = u1 .- u0
        calculate_residuals(ũ, u0, u1, α, ρ, ODE_DEFAULT_NORM, t)
        calculate_residuals!(out, ũ, u0, u1, α, ρ, ODE_DEFAULT_NORM, t)

        # Precompile scalar versions (commonly used)
        ODE_DEFAULT_NORM(1.0, t)
        NAN_CHECK(1.0)
        INFINITE_OR_GIANT(1.0)
        ODE_DEFAULT_UNSTABLE_CHECK(dt, 1.0, nothing, t)
    end
end
