using Test
using SparseArrays
using OrdinaryDiffEqNonlinearSolve: find_algebraic_vars_eqs

@testset "Sparse Algebraic Detection Performance" begin
    # Test 1: Correctness - results should match between sparse and dense methods
    @testset "Correctness" begin
        # Small mixed matrix
        M1 = sparse([1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 2.0])
        vars1, eqs1 = find_algebraic_vars_eqs(M1)
        @test vars1 == [false, true, false]  # col 2 is zero (algebraic)
        @test eqs1 == [false, true, false]   # row 2 is zero (algebraic)

        # Diagonal matrix (no algebraic)
        M2 = sparse(1:10, 1:10, ones(10))
        vars2, eqs2 = find_algebraic_vars_eqs(M2)
        @test all(.!vars2)  # no algebraic vars
        @test all(.!eqs2)   # no algebraic eqs

        # Zero matrix (all algebraic)
        M3 = spzeros(5, 5)
        vars3, eqs3 = find_algebraic_vars_eqs(M3)
        @test all(vars3)  # all algebraic vars
        @test all(eqs3)   # all algebraic eqs

        # Matrix with explicit zeros in sparsity pattern
        M4 = sparse([1, 2], [1, 2], [1.0, 0.0], 3, 3)  # (2,2) is stored but zero
        vars4, eqs4 = find_algebraic_vars_eqs(M4)
        @test vars4 == [false, true, true]   # cols 2,3 are effectively zero
        @test eqs4 == [false, true, true]    # rows 2,3 are effectively zero
    end

    # Test 2: Complexity verification
    # For sparse matrices, time should scale with nnz, not n²
    # We verify this by checking that processing a large sparse diagonal matrix
    # is much faster than it would be with O(n²) complexity
    @testset "Complexity O(nnz) not O(n²)" begin
        # Create diagonal matrices of increasing size
        # For a diagonal matrix: nnz = n, so O(nnz) = O(n)
        # If algorithm were O(n²), doubling n would 4x the time
        # With O(nnz) = O(n), doubling n should ~2x the time

        n_small = 10_000
        n_large = 40_000

        M_small = sparse(1:n_small, 1:n_small, ones(n_small))
        M_large = sparse(1:n_large, 1:n_large, ones(n_large))

        # Warmup
        find_algebraic_vars_eqs(M_small)
        find_algebraic_vars_eqs(M_large)

        # Measure
        t_small = @elapsed for _ in 1:10
            find_algebraic_vars_eqs(M_small)
        end

        t_large = @elapsed for _ in 1:10
            find_algebraic_vars_eqs(M_large)
        end

        ratio = t_large / t_small
        size_ratio = n_large / n_small  # 4x size increase

        # With O(n²), ratio would be ~16
        # With O(nnz) = O(n), ratio should be ~4
        # Allow some tolerance for measurement noise
        @test ratio < size_ratio * 2  # Should be closer to linear than quadratic

        # Absolute sanity check: processing 40k diagonal shouldn't take > 100ms
        t_single = @elapsed find_algebraic_vars_eqs(M_large)
        @test t_single < 0.1  # Should be < 100ms, typically < 1ms
    end

    # Test 3: Verify fallback for dense matrices still works
    @testset "Dense matrix fallback" begin
        M_dense = [1.0 0.0; 0.0 0.0]
        vars, eqs = find_algebraic_vars_eqs(M_dense)
        @test vars == [false, true]
        @test eqs == [false, true]
    end
end
