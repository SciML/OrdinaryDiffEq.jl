@testset "Implicit Taylor stage accuracy" begin
    ODTI = OrdinaryDiffEqTaylorSeriesImplicit
    scaled_norm(v, scale) = sqrt(sum(abs2, v ./ scale) / length(v))

    # Follow the branch emanating from h=0, independently of the numerical
    # output. Agreement after halving the continuation spacing checks that
    # these sampled reference solves have not jumped to another stage root.
    function continued_root(root, u0, h, alg, subdivisions)
        T = alg.μ isa Complex ? Complex{BigFloat} : BigFloat
        x = T.(u0)
        for j in 1:subdivisions
            x = root(u0, x, BigFloat(h) * j / subdivisions, alg)
        end
        return x
    end

    @testset "Stage accuracy under changes of state units" begin
        # Independent scalar Taylor coefficients for y' = lambda*y^2/s:
        # c_k(y) = (lambda/s)^k*y^(k+1). These also give a reference for
        # complex mu-Taylor stages without differentiating a generated kernel.
        function riccati_root(u0, guess, h, alg, scales)
            T = eltype(guess) <: Complex ? Complex{BigFloat} : BigFloat
            x, v = T.(guess), T.(u0)
            P, Q = ODTI._unwrap_val(alg.order), ODTI._unwrap_val(alg.order_q)
            rates = T.((1, -1)) ./ BigFloat.(scales)
            rhs = zero(x)
            if ODTI.is_taylor_pade(alg)
                a, b = ODTI.pade_increment_coefficients(Val(P), Val(Q))
                d = BigFloat(h)
                for k in 1:P
                    rhs .+= T(a[k + 1]) * d^k .* rates.^k .* v.^(k + 1)
                end
            else
                b = ODTI.taylor_increment_coefficients(Val(P))
                d = -T(alg.μ) * BigFloat(h)
                Q = P
            end
            for _ in 1:30
                residual = x - v - rhs
                derivative = ones(T, length(x))
                for k in 1:Q
                    w = T(b[k + 1]) * d^k .* rates.^k
                    residual .+= w .* x.^(k + 1)
                    derivative .+= (k + 1) .* w .* x.^k
                end
                correction = residual ./ derivative
                x -= correction
                maximum(abs.(correction ./ scales)) < big"1e-50" && return x
            end
            error("Independent Riccati stage solve did not converge")
        end
        riccati!(du, u, scales, t) = begin
            du[1] = u[1]^2 / scales[1]
            du[2] = -u[2]^2 / scales[2]
            nothing
        end
        methods = (TaylorGaussA2, TaylorGaussA4, TaylorRadauLA3,
            TaylorLobattoLA4, TaylorRadauLA5, ImplicitTaylor1A2, ImplicitTaylor2A4)
        for make in methods
            normalized = []
            decisions = []
            for scales in ([1.0, 1.0], [1e-9, 1e6])
                u0 = scales .* [0.5, 1.0]
                problem = ODEProblem(riccati!, u0, (0.0, 0.5), scales)
                i = init(problem, make(); dt = 0.1, abstol = 1e-8 .* scales, reltol = 1e-5)
                @test ODTI.implicit_taylor_step!(i, i.cache, 0.0, 0.1, u0, i.u)
                root(args...) = riccati_root(args..., scales)
                r32 = continued_root(root, u0, 0.1, i.alg, 32)
                r64 = continued_root(root, u0, 0.1, i.alg, 64)
                @test maximum(abs.((r32 - r64) ./ scales)) < big"1e-35"
                scale = 1e-8 .* scales .+ 1e-5 .* max.(abs.(u0), abs.(i.cache.uintermediate))
                @test scaled_norm(i.cache.uintermediate - r64, scale) < i.cache.κ
                push!(normalized, i.cache.uintermediate ./ scales)
                push!(decisions, OrdinaryDiffEqCore.get_EEst(i) <= 1)
            end
            @test normalized[1] ≈ normalized[2] rtol = 1e-10 atol = 1e-13
            @test decisions[1] == decisions[2]
        end
        for scales in ([1.0, 1.0], [1e-9, 1e6]), make in methods
            u0 = scales .* [0.5, 1.0]
            problem = ODEProblem(riccati!, u0, (0.0, 0.5), scales)
            abstol, reltol = 1e-8 .* scales, 1e-5
            i = init(problem, make(); dt = 0.1, abstol, reltol)
            for _ in 1:1000
                i.t >= 0.5 && break
                old, t0 = copy(i.u), i.t
                step!(i)
                @test i.t > t0
                c = i.cache
                reference = riccati_root(old, c.uintermediate, i.dt, i.alg, scales)
                scale = abstol .+ reltol .* max.(abs.(old), abs.(c.uintermediate))
                error = scaled_norm(c.uintermediate - reference, scale)
                @test error < c.κ
            end
            @test i.t == 0.5
        end
    end

    @testset "Quadratic stiff kinetics stage roots" begin
        # Reconstruct the jet from the original quadratic RHS, independently
        # of Symbolics and the solver's coefficient-Jacobian assembly:
        # c1=f, c2=J*f/2, c3=(J*J*f+H[f,f])/6.
        problem = ODEProblemLibrary.prob_ode_pollution
        rhs! = SciMLBase.unwrapped_f(problem.f)
        rate(u) = (du = zero(u); rhs!(du, u, problem.p, zero(eltype(u))); du)
        n = length(problem.u0)
        L = ForwardDiff.jacobian(rate, zeros(BigFloat, n))
        H = cat([ForwardDiff.hessian(u -> rate(u)[i], zeros(BigFloat, n))
                 for i in 1:n]...; dims = 3)
        nonzeros = [(i, j, k, H[j, k, i]) for i in 1:n, j in 1:n, k in 1:n
                    if !iszero(H[j, k, i])]
        function hessian_action(v)
            A = zeros(eltype(v), n, n)
            for (i, j, k, weight) in nonzeros
                A[i, j] += weight * v[k]
            end
            return A
        end
        function quadratic_jet(v)
            f = rate(v)
            J = L + hessian_action(v)
            g, Hf = J * f, hessian_action(f)
            coefficients = (f, g / 2, (J * g + Hf * f) / 6)
            derivatives = (J, (Hf + J * J) / 2,
                (hessian_action(g) + J * Hf + J * J * J + 2Hf * J) / 6)
            return coefficients, derivatives
        end
        function quadratic_root(u0, guess, h, alg)
            P, Q = ODTI._unwrap_val(alg.order), ODTI._unwrap_val(alg.order_q)
            a, b = ODTI.pade_increment_coefficients(Val(P), Val(Q))
            v, x, d = BigFloat.(u0), BigFloat.(guess), BigFloat(h)
            coefficients, _ = quadratic_jet(v)
            rhs = sum(BigFloat(a[k + 1]) * d^k * coefficients[k] for k in 1:P)
            for _ in 1:30
                coefficients, derivatives = quadratic_jet(x)
                residual = x - v - rhs
                M = Matrix{BigFloat}(I, n, n)
                for k in 1:Q
                    weight = BigFloat(b[k + 1]) * d^k
                    residual += weight * coefficients[k]
                    M += weight * derivatives[k]
                end
                correction = M \ residual
                x -= correction
                maximum(abs.(correction)) < big"1e-40" && return x
            end
            error("Independent quadratic stage solve did not converge")
        end
        for make in (TaylorLobattoLA4, TaylorRadauLA5), reltol in (1e-3, 1e-5)
            abstol = reltol / 1000
            i = init(problem, make(); abstol, reltol, save_everystep = false)
            worst = 0.0
            checkpoints = (0.0, 0.1, 1.0, 5.0, 14.0, 40.0)
            checkpoint = 1
            for _ in 1:2000
                i.t >= problem.tspan[2] && break
                old, t0 = copy(i.u), i.t
                step!(i)
                @test i.t > t0
                reference = quadratic_root(old, i.u, i.dt, i.alg)
                scale = abstol .+ reltol .* max.(abs.(old), abs.(i.u))
                error = scaled_norm(i.u - reference, scale)
                worst = max(worst, error)
                @test error < i.cache.κ
                if checkpoint <= length(checkpoints) && t0 >= checkpoints[checkpoint]
                    r32 = continued_root(quadratic_root, old, i.dt, i.alg, 32)
                    r64 = continued_root(quadratic_root, old, i.dt, i.alg, 64)
                    @test maximum(abs.(r32 - r64)) < big"1e-30"
                    @test scaled_norm(i.u - r64, scale) < i.cache.κ
                    checkpoint += 1
                end
            end
            @test i.t == problem.tspan[2]
            @test worst < i.cache.κ
            @test checkpoint == length(checkpoints) + 1
        end
    end
end
