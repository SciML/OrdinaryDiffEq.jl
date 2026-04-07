"""
    deduce_Butcher_tableau(erk, T = Float64)

Deduce and return the Butcher coefficients `A, b, c` by solving some
specific ordinary differential equations using the explicit Runge-Kutta
method `erk`. The type `T` will be used for computations and is the
`eltype` of `A`, `b`, and `c`.
"""
function deduce_Butcher_tableau(erk, T = Float64)
    # get the number of stages s
    step_counter_s = 0
    function f_for_s(du, u, p, t)
        step_counter_s += 1
        du .= zero(eltype(du))
        nothing
    end

    tspan = (zero(T), one(T))
    u0 = [zero(T)]
    ode = ODEProblem(f_for_s, u0, tspan)
    step_counter_s = 0
    integrator = init(ode, erk, adaptive = false, dt = 1.0)
    step!(integrator)
    s = step_counter_s

    # get the nodes c[i]
    step_counter_c = 0
    c = zeros(T, s)

    function f_for_c(du, u, p, t)
        step_counter_c += 1
        if step_counter_c <= s
            c[step_counter_c] = t
        end
        du .= zero(eltype(du))
        nothing
    end

    tspan = (zero(T), one(T))
    u0 = [zero(T)]
    ode = ODEProblem(f_for_c, u0, tspan)
    step_counter_c = 0
    integrator = init(ode, erk, adaptive = false, dt = 1.0)
    step!(integrator)

    # get the weights b[i]
    step_counter_b = 0
    b = zeros(T, s)
    function f_for_b(du, u, p, t)
        step_counter_b += 1
        for idx in 1:length(du)
            du[idx] = step_counter_b == idx
        end
        nothing
    end

    tspan = (zero(T), one(T))
    u0 = zeros(T, s)
    ode = ODEProblem(f_for_b, u0, tspan)
    step_counter_b = 0
    integrator = init(ode, erk, adaptive = false, dt = 1.0)
    step!(integrator)
    for i in 1:s
        b[i] = integrator.u[i]
    end

    # get the coefficients A[i,j]
    step_counter_A = 0
    A = zeros(T, s, s)
    function f_for_A(du, u, p, t)
        step_counter_A += 1
        for idx in 1:length(du)
            A[step_counter_A, idx] = u[idx]
            du[idx] = step_counter_A == idx
        end
        nothing
    end

    tspan = (zero(T), one(T))
    u0 = zeros(T, s)
    ode = ODEProblem(f_for_A, u0, tspan)
    step_counter_A = 0
    integrator = init(ode, erk, adaptive = false, dt = 1.0)
    step!(integrator)

    A, b, c
end
