using Test
using SparseArrays
using OrdinaryDiffEqCore: _isdiag

@testset "Sparse isdiag Performance" begin
    # Test 1: Correctness - _isdiag should correctly identify diagonal matrices
    @testset "Correctness" begin
        # Diagonal sparse matrix
        M1 = sparse(1:5, 1:5, ones(5))
        @test _isdiag(M1) == true

        # Non-diagonal sparse matrix (has off-diagonal elements)
        M2 = sparse([1, 1, 2], [1, 2, 2], [1.0, 2.0, 3.0], 3, 3)
        @test _isdiag(M2) == false

        # Empty sparse matrix (is diagonal by definition)
        M3 = spzeros(5, 5)
        @test _isdiag(M3) == true

        # Non-square matrix (not diagonal)
        M4 = sparse([1, 2], [1, 2], [1.0, 2.0], 3, 4)
        @test _isdiag(M4) == false

        # Sparse matrix with explicit zeros stored on diagonal
        M5 = sparse([1, 2, 3], [1, 2, 3], [1.0, 0.0, 1.0], 3, 3)
        @test _isdiag(M5) == true

        # Large diagonal matrix with some zeros
        n = 100
        vals = vcat(ones(n ÷ 2), zeros(n - n ÷ 2))
        M6 = spdiagm(0 => vals)
        @test _isdiag(M6) == true
    end

    # Test 2: Complexity verification
    # For sparse matrices, time should scale with nnz, not n^2
    @testset "Complexity O(nnz) not O(n^2)" begin
        n_small = 10_000
        n_large = 40_000

        M_small = sparse(1:n_small, 1:n_small, ones(n_small))
        M_large = sparse(1:n_large, 1:n_large, ones(n_large))

        # Warmup
        _isdiag(M_small)
        _isdiag(M_large)

        # Measure
        t_small = @elapsed for _ in 1:10
            _isdiag(M_small)
        end

        t_large = @elapsed for _ in 1:10
            _isdiag(M_large)
        end

        ratio = t_large / t_small
        size_ratio = n_large / n_small  # 4x size increase

        # With O(n^2), ratio would be ~16
        # With O(nnz) = O(n), ratio should be ~4
        # Allow some tolerance for measurement noise
        @test ratio < size_ratio * 2  # Should be closer to linear than quadratic

        # Absolute sanity check: processing 40k diagonal shouldn't take > 100ms
        t_single = @elapsed _isdiag(M_large)
        @test t_single < 0.1  # Should be < 100ms, typically < 1ms
    end

    # Test 3: Verify fallback for dense matrices still works
    @testset "Dense matrix fallback" begin
        M_dense = [1.0 0.0; 0.0 2.0]
        @test _isdiag(M_dense) == true

        M_dense_nondiag = [1.0 1.0; 0.0 2.0]
        @test _isdiag(M_dense_nondiag) == false
    end
end
