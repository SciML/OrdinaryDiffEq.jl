using StochasticDiffEq, DiffEqNoiseProcess, SparseArrays, LinearAlgebra,
    AllocCheck
using DiffEqBase: @..

@testset "EulerHeun sparse noise: no per-step alloc" begin

    # Simple linear drift
    f!(du, u, p, t) = (@.. du = 0.999 * u)

    # 2×2 identical-column block structure; g! only writes nzval of an existing sparsity pattern
    function sparse_proto(N)
        I = Vector{Int}(undef, 4N)
        J = similar(I)
        V = ones(Float64, 4N)
        @inbounds for i in 1:N
            I[4i - 3] = 2i - 1
            J[4i - 3] = 2i - 1

            I[4i - 2] = 2i
            J[4i - 2] = 2i - 1

            I[4i - 1] = 2i - 1
            J[4i - 1] = 2i

            I[4i] = 2i
            J[4i] = 2i
        end
        sparse(I, J, V, 2N, 2N)
    end

    @inline function ensure_pattern!(G::SparseMatrixCSC{T, Int}, N) where {T}
        s = one(T)
        @inbounds for i in 1:N
            G[2i - 1, 2i - 1] = s
            G[2i, 2i - 1] = s
            G[2i - 1, 2i] = s
            G[2i, 2i] = s
        end
        return nothing
    end

    # Dense g!
    function g!(G::StridedMatrix{T}, u, p, t) where {T}
        c012 = T(0.12)
        c18 = T(1.8)
        @inbounds for i in 1:(p.N)
            off = 2i - 1
            G[off, off] = c012 * u[2i - 1]
            G[off + 1, off] = c18 * u[2i]

            G[off, off + 1] = c012 * u[2i - 1]
            G[off + 1, off + 1] = c18 * u[2i]
        end
        return nothing
    end

    # Sparse g!
    function g!(G::SparseMatrixCSC{T}, u, p, t) where {T}
        c012 = T(0.12)
        c18 = T(1.8)
        @inbounds for i in 1:(p.N)
            off = G.colptr[2i - 1]
            G.nzval[off] = c012 * u[2i - 1]
            G.nzval[off + 1] = c18 * u[2i]

            off = G.colptr[2i]
            G.nzval[off] = c012 * u[2i - 1]
            G.nzval[off + 1] = c18 * u[2i]
        end
        return nothing
    end

    function make_integrator_sparse(N)
        A = sparse_proto(N)
        p = (; N)
        W = SimpleWienerProcess!(0.0, zeros(2N); save_everystep = false)
        prob = SDEProblem(
            f!, g!, ones(2N), (0.0, 1.0), p; noise_rate_prototype = A, noise = W
        )
        integ = init(prob, EulerHeun(); dt = 0.01, adaptive = false, save_on = false)

        cache = integ.cache
        allocs = AllocCheck.check_allocs(
            StochasticDiffEq.perform_step!, (typeof(integ), typeof(cache))
        )
        @test isempty(allocs)
    end

    function make_integrator_dense(N)
        A = zeros(2N, 2N)
        p = (; N)
        W = SimpleWienerProcess!(0.0, zeros(2N); save_everystep = false)
        prob = SDEProblem(
            f!, g!, ones(2N), (0.0, 1.0), p; noise_rate_prototype = A, noise = W
        )
        integ = init(prob, EulerHeun(); dt = 0.01, adaptive = false, save_on = false)

        cache = integ.cache

        # Dense+BLAS: assert with `@allocated == 0` (not `check_allocs`).
        # `check_allocs` flags throw-only branches in LinearAlgebra’s generic gemv!/mul!,
        # while the BLAS hot path is allocation-free at runtime. We keep `check_allocs`
        # for the sparse/non-diagonal tests.
        StochasticDiffEq.perform_step!(integ, cache) # warm-up
        @test @allocated(StochasticDiffEq.perform_step!(integ, cache)) == 0
    end

    make_integrator_dense(16)
    make_integrator_sparse(16)
end
