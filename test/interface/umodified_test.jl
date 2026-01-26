using OrdinaryDiffEq, Test, LinearAlgebra

T = 1000.0;
Ttr = 0.0;
d0 = 1.0e-9;
threshold = 10^4 * d0;
dt = 0.1;
diff_eq_kwargs = Dict();

@inbounds function eom_lorenz!(du, u, p, t)
    σ = 10.0
    ρ = 28.0
    β = 8 / 3
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - β * u[3]
end

prob = ODEProblem(eom_lorenz!, [0.0, 10.0, 0], (zero(T), T))
integ1 = init(prob, Tsit5(); diff_eq_kwargs...)
@. prob.u0 = prob.u0 + d0 / sqrt(3)
integ2 = init(prob, Tsit5(); diff_eq_kwargs...)
integ1.opts.advance_to_tstop = true
integ2.opts.advance_to_tstop = true

dist = d0
λ = zero(eltype(integ1.u))
i = 0
tvector = dt:dt:T

for τ in tvector
    # evolve until rescaling:
    push!(integ1.opts.tstops, τ)
    step!(integ1)
    push!(integ2.opts.tstops, τ)
    step!(integ2)
    global dist = norm(integ1.u .- integ2.u)
    # Rescale:
    if dist ≥ threshold
        # Rescale and reset everything:
        integ2.u .= integ1.u #.+ (integ2.u .- integ1.u)./a
        u_modified!(integ2, true)
        set_proposed_dt!(integ2, integ1)
        if hasproperty(integ1, :controller_cache)
            reinit!(integ1, integ1.controller_cache)
            reinit!(integ2, integ2.controller_cache)
        end
        break
    end
end

τ = tvector[end]
push!(integ1.opts.tstops, τ);
step!(integ1);
@test integ1.t == τ
push!(integ2.opts.tstops, τ);
step!(integ2);
@test integ2.t == τ
@test dist = norm(integ1.u .- integ2.u) == 0
