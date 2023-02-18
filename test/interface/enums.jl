using OrdinaryDiffEq, PoissonRandom

function rate_to_proportion(r::Float64, t::Float64)
    1 - exp(-r * t)
end;

function sir_abm!(du, u, p, t)
    (β, c, γ, δt) = p
    N = length(u)
    # Initialize du to u
    for i in 1:N
        du[i] = u[i]
    end
    for i in 1:N # loop through agents
        # If recovered
        if u[i] == Recovered
            continue
            # If susceptible
        elseif u[i] == Susceptible
            ncontacts = pois_rand(c * δt)
            while ncontacts > 0
                j = rand(1:N)
                if j == i
                    continue
                end
                a = u[j]
                if a == Infected && rand() < β
                    du[i] = Infected
                    break
                end
                ncontacts -= 1
            end
            # If infected
        else
            u[i] == Infected
            if rand() < γ
                du[i] = Recovered
            end
        end
    end
    nothing
end

δt = 0.1
nsteps = 400
tf = nsteps * δt
tspan = (0.0, nsteps)
t = 0:δt:tf;
β = 0.05
c = 10.0
γ = rate_to_proportion(0.25, δt)
p = [β, c, γ, δt]

@enum InfectionStatus Susceptible Infected Recovered

N = 1000
I0 = 10
u0 = Array{InfectionStatus}(undef, N);
for i in 1:N
    if i <= I0
        s = Infected
    else
        s = Susceptible
    end
    u0[i] = s
end

prob_map = DiscreteProblem(sir_abm!, u0, tspan, p)
sol_map = solve(prob_map, FunctionMap())
