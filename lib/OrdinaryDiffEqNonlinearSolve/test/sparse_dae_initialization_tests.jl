using Test
using SparseArrays
using LinearAlgebra
using OrdinaryDiffEqSDIRK
using OrdinaryDiffEqNonlinearSolve: find_algebraic_vars_eqs, algebraic_jacobian
using SciMLBase: ReturnCode

@testset "Sparse jac_prototype in BrownFullBasicInit" begin
    # Verify algebraic_jacobian preserves sparsity
    @testset "algebraic_jacobian preserves SparseMatrixCSC" begin
        N = 20
        # Build a sparse jac_prototype for N differential + N algebraic
        I_idx = Int[]
        J_idx = Int[]
        for i in 1:N
            push!(I_idx, i); push!(J_idx, i)       # diff: du_i/du_i
            push!(I_idx, i); push!(J_idx, N + i)    # diff: du_i/dalg_i
            push!(I_idx, N + i); push!(J_idx, i)    # alg: dalg_i/du_i
            push!(I_idx, N + i); push!(J_idx, N + i) # alg: dalg_i/dalg_i
        end
        jac_proto = sparse(I_idx, J_idx, ones(length(I_idx)), 2N, 2N)

        # Diagonal mass matrix
        M_diag = Diagonal(vcat(ones(N), zeros(N)))
        alg_vars_d, alg_eqs_d = find_algebraic_vars_eqs(M_diag)
        J_d = algebraic_jacobian(jac_proto, alg_eqs_d, alg_vars_d)
        @test J_d isa SparseMatrixCSC
        @test size(J_d) == (N, N)
        @test nnz(J_d) == N  # only diagonal coupling in algebraic block

        # Sparse mass matrix
        M_sparse = sparse(M_diag)
        alg_vars_s, alg_eqs_s = find_algebraic_vars_eqs(M_sparse)
        J_s = algebraic_jacobian(jac_proto, alg_eqs_s, alg_vars_s)
        @test J_s isa SparseMatrixCSC
        @test J_s == J_d

        # nothing jac_prototype returns nothing
        @test algebraic_jacobian(nothing, alg_eqs_d, alg_vars_d) === nothing
    end

    # Small system: 1 differential + 1 algebraic (in-place)
    @testset "Small in-place ODEProblem with sparse jac_prototype" begin
        function f_small!(du, u, p, t)
            du[1] = -u[1] + u[2]
            du[2] = u[1] - u[2]^3
            nothing
        end
        M = [1.0 0.0; 0.0 0.0]
        jac_proto = sparse([1.0 1.0; 1.0 1.0])

        u0 = [1.0, 0.5]  # inconsistent (u1 - u2^3 ≠ 0)
        f_ode = ODEFunction(f_small!; mass_matrix = M, jac_prototype = jac_proto)
        prob = ODEProblem(f_ode, u0, (0.0, 1.0))

        sol = solve(prob, Trapezoid())
        @test sol.retcode == ReturnCode.Success
        # Check constraint is satisfied after initialization
        @test abs(sol.u[1][1] - sol.u[1][2]^3) < 1.0e-6
    end

    # Out-of-place: test algebraic_jacobian produces correct sparse sub-matrix
    @testset "Out-of-place algebraic_jacobian slicing" begin
        M = [1.0 0.0; 0.0 0.0]
        jac_proto = sparse([1.0 1.0; 1.0 1.0])
        alg_vars, alg_eqs = find_algebraic_vars_eqs(M)
        J = algebraic_jacobian(jac_proto, alg_eqs, alg_vars)
        @test J isa SparseMatrixCSC
        @test size(J) == (1, 1)
        @test J[1, 1] == 1.0
    end

    # ROBER problem with sparse jac_prototype
    @testset "ROBER with sparse jac_prototype" begin
        function rober!(du, u, p, t)
            y₁, y₂, y₃ = u
            k₁, k₂, k₃ = p
            du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
            du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
            du[3] = y₁ + y₂ + y₃ - 1
            nothing
        end
        M = [1.0 0 0; 0 1.0 0; 0 0 0]
        # Full sparsity for ROBER (all entries nonzero)
        jac_proto = sparse(ones(3, 3))
        p = (0.04, 3.0e7, 1.0e4)

        # Inconsistent initial condition
        u0 = [1.0, 0.0, 0.2]
        f_ode = ODEFunction(rober!; mass_matrix = M, jac_prototype = jac_proto)
        prob = ODEProblem(f_ode, u0, (0.0, 1.0e5), p)

        sol = solve(prob, Trapezoid())
        @test sol.retcode == ReturnCode.Success
        @test sum(sol.u[1]) ≈ 1
        @test sol.u[1] ≈ [1.0, 0.0, 0.0]
    end

    # Medium-sized sparse system (N diff + N alg with nearest-neighbor coupling)
    @testset "Medium sparse system (N=$N)" for N in [50, 200]
        function f_medium!(du, u, p, t)
            n = div(length(u), 2)
            for i in 1:n
                prev = i > 1 ? u[i - 1] : zero(eltype(u))
                next = i < n ? u[i + 1] : zero(eltype(u))
                du[i] = -2 * u[i] + prev + next + u[n + i]
            end
            for i in 1:n
                du[n + i] = u[i] - u[n + i]^2
            end
            nothing
        end

        M = Diagonal(vcat(ones(N), zeros(N)))

        # Build sparse jac_prototype: tridiagonal for diff + diagonal coupling
        I_idx = Int[]
        J_idx = Int[]
        for i in 1:N
            push!(I_idx, i); push!(J_idx, i)
            if i > 1
                push!(I_idx, i); push!(J_idx, i - 1)
            end
            if i < N
                push!(I_idx, i); push!(J_idx, i + 1)
            end
            push!(I_idx, i); push!(J_idx, N + i)
            push!(I_idx, N + i); push!(J_idx, i)
            push!(I_idx, N + i); push!(J_idx, N + i)
        end
        jac_proto = sparse(I_idx, J_idx, ones(length(I_idx)), 2N, 2N)

        u0 = vcat(ones(N), 0.5 * ones(N))  # inconsistent
        f_ode = ODEFunction(f_medium!; mass_matrix = M, jac_prototype = jac_proto)
        prob = ODEProblem(f_ode, u0, (0.0, 1.0))

        sol = solve(prob, Trapezoid(); save_everystep = false)
        @test sol.retcode == ReturnCode.Success
        # Algebraic constraint u[i] = u[N+i]^2 should be satisfied at t=0
        max_violation = maximum(abs.(sol.u[1][1:N] .- sol.u[1][(N + 1):(2N)] .^ 2))
        @test max_violation < 1.0e-5
    end

    # Verify sparse mass matrix works the same way
    @testset "Sparse mass matrix with sparse jac_prototype" begin
        N = 30
        function f_sp!(du, u, p, t)
            n = div(length(u), 2)
            for i in 1:n
                du[i] = -u[i] + u[n + i]
            end
            for i in 1:n
                du[n + i] = u[i] - u[n + i]^2
            end
            nothing
        end

        M_dense = Diagonal(vcat(ones(N), zeros(N)))
        M_sparse = sparse(M_dense)

        I_idx = Int[]
        J_idx = Int[]
        for i in 1:N
            push!(I_idx, i); push!(J_idx, i)
            push!(I_idx, i); push!(J_idx, N + i)
            push!(I_idx, N + i); push!(J_idx, i)
            push!(I_idx, N + i); push!(J_idx, N + i)
        end
        jac_proto = sparse(I_idx, J_idx, ones(length(I_idx)), 2N, 2N)

        u0 = vcat(ones(N), 0.5 * ones(N))

        f_dense = ODEFunction(f_sp!; mass_matrix = M_dense, jac_prototype = jac_proto)
        f_sparse = ODEFunction(f_sp!; mass_matrix = M_sparse, jac_prototype = jac_proto)

        prob_dense = ODEProblem(f_dense, u0, (0.0, 1.0))
        prob_sparse = ODEProblem(f_sparse, u0, (0.0, 1.0))

        sol_dense = solve(prob_dense, Trapezoid(); save_everystep = false)
        sol_sparse = solve(prob_sparse, Trapezoid(); save_everystep = false)

        @test sol_dense.retcode == ReturnCode.Success
        @test sol_sparse.retcode == ReturnCode.Success
        @test sol_dense.u[1] ≈ sol_sparse.u[1]
    end

    # Large-scale test inspired by issue DifferentialEquations.jl#1107
    # (2D grid, two coupled fields)
    @testset "Large sparse DAE (issue #1107 pattern)" begin
        # Nx x Ny grid with 2 fields: P (differential) and phi (algebraic)
        Nx, Ny = 32, 16
        N = Nx * Ny  # unknowns per field
        ntotal = 2N  # total unknowns

        # Linear index helper
        idx(i, j) = (j - 1) * Nx + i
        idx_P(i, j) = idx(i, j)
        idx_phi(i, j) = N + idx(i, j)

        function f_2d!(du, u, p, t)
            for j in 1:Ny, i in 1:Nx
                # Differential: dP/dt = laplacian(P) + phi
                ip = i < Nx ? i + 1 : 1
                im = i > 1 ? i - 1 : Nx
                jp = j < Ny ? j + 1 : 1
                jm = j > 1 ? j - 1 : Ny
                lap = u[idx_P(ip, j)] + u[idx_P(im, j)] +
                    u[idx_P(i, jp)] + u[idx_P(i, jm)] - 4 * u[idx_P(i, j)]
                du[idx_P(i, j)] = lap + u[idx_phi(i, j)]

                # Algebraic: 0 = P - phi^3
                du[idx_phi(i, j)] = u[idx_P(i, j)] - u[idx_phi(i, j)]^3
            end
            nothing
        end

        # Mass matrix: I for P equations, 0 for phi equations
        M = Diagonal(vcat(ones(N), zeros(N)))

        # Build sparse jac_prototype
        I_idx = Int[]
        J_idx = Int[]
        for j in 1:Ny, i in 1:Nx
            k = idx(i, j)
            # dP equation depends on 5-point stencil of P + phi(i,j)
            push!(I_idx, k); push!(J_idx, k)  # P(i,j)
            ip = i < Nx ? i + 1 : 1
            im = i > 1 ? i - 1 : Nx
            jp = j < Ny ? j + 1 : 1
            jm = j > 1 ? j - 1 : Ny
            push!(I_idx, k); push!(J_idx, idx(ip, j))
            push!(I_idx, k); push!(J_idx, idx(im, j))
            push!(I_idx, k); push!(J_idx, idx(i, jp))
            push!(I_idx, k); push!(J_idx, idx(i, jm))
            push!(I_idx, k); push!(J_idx, N + k)  # phi(i,j)

            # phi equation depends on P(i,j) and phi(i,j)
            push!(I_idx, N + k); push!(J_idx, k)
            push!(I_idx, N + k); push!(J_idx, N + k)
        end
        jac_proto = sparse(I_idx, J_idx, ones(length(I_idx)), ntotal, ntotal)

        # Verify sparsity structure
        @test nnz(jac_proto) < ntotal * ntotal * 0.01  # very sparse

        # Inconsistent initial condition
        u0 = vcat(ones(N), 0.5 * ones(N))
        f_ode = ODEFunction(f_2d!; mass_matrix = M, jac_prototype = jac_proto)
        prob = ODEProblem(f_ode, u0, (0.0, 0.1))

        # Verify the algebraic sub-jacobian is sparse and correct size
        alg_vars, alg_eqs = find_algebraic_vars_eqs(M)
        J_sub = algebraic_jacobian(jac_proto, alg_eqs, alg_vars)
        @test J_sub isa SparseMatrixCSC
        @test size(J_sub) == (N, N)
        @test nnz(J_sub) == N  # diagonal only for algebraic block

        sol = solve(prob, Trapezoid(); save_everystep = false)
        @test sol.retcode == ReturnCode.Success
        # Verify algebraic constraint P = phi^3 after initialization
        max_violation = maximum(
            abs.(sol.u[1][1:N] .- sol.u[1][(N + 1):(2N)] .^ 3)
        )
        @test max_violation < 1.0e-4
    end
end
