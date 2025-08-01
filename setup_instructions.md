# OrdinaryDiffEq.jl Precompilation Setup Instructions

## Files Created

1. **`precompile_workload.jl`** - Comprehensive precompilation script
2. **`CI_with_precompile.yml`** - GitHub Actions workflow with precompile job
3. **Setup instructions** (this file)

## Integration Steps

### 1. Add Precompilation Script to Repository

Copy `precompile_workload.jl` to the root of the OrdinaryDiffEq.jl repository:

```bash
cp precompile_workload.jl /path/to/OrdinaryDiffEq.jl/
```

### 2. Add PrecompileTools Dependency

Add PrecompileTools to the main Project.toml:

```toml
[deps]
# ... existing dependencies ...
PrecompileTools = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
```

### 3. Update GitHub Actions Workflow

Replace or modify the existing `.github/workflows/CI.yml` with the new workflow, or integrate the precompile job pattern:

```bash
cp CI_with_precompile.yml /path/to/OrdinaryDiffEq.jl/.github/workflows/
```

### 4. Key Benefits

- **Shared Precompilation**: One precompile job builds cache, all test jobs reuse it
- **Faster CI**: Test jobs skip expensive first-compilation of core algorithms  
- **Comprehensive Coverage**: Exercises Tsit5, FBDF, Rodas5P, Vern9 plus common patterns
- **Production-Ready**: Same precompilation benefits users in real applications

### 5. Workflow Structure

```
precompile (builds cache)
├── test-interface (depends on precompile)
├── test-integrators (depends on precompile)  
├── test-algorithms (matrix, depends on precompile)
├── test-subpackages (matrix, depends on precompile)
├── test-performance (depends on precompile)
├── test-multithreading (depends on precompile)
└── coverage (depends on test jobs, runs with coverage=true)
```

### 6. Cache Strategy

- Uses `julia-actions/cache@v2` with consistent `cache-name: 'precompile-cache'`
- Precompile job populates cache, all dependent jobs reuse it
- Cross-platform jobs get their own cache keys to avoid conflicts

### 7. Customization Options

#### Modify Precompilation Scope
Edit `precompile_workload.jl` to add/remove solver algorithms or problem types.

#### Adjust Test Matrix
Modify the `strategy.matrix` sections to match your current test groups.

#### Add More Solvers
Extend the precompilation script with additional algorithms:

```julia
println("  - Precompiling RadauIIA5 (implicit Runge-Kutta)...")
solve(prob_vdp, RadauIIA5(), reltol=1e-8, abstol=1e-10)

println("  - Precompiling TRBDF2 (trapezoidal BDF)...")  
solve(prob_stiff, TRBDF2(), reltol=1e-6, abstol=1e-8)
```

## Performance Expectations

Based on the SciML compile time analysis, this setup should:

- Reduce first-solve time from ~22s to ~3s for common algorithms
- Speed up CI by 2x-4x for test jobs (excluding the precompile job time)
- Improve developer experience with faster local test runs
- Benefit production users with pre-warmed solver compilation

## Monitoring Success

Check CI logs for:
- Precompile job completes successfully
- Test jobs show faster startup times  
- Cache hit rates in julia-actions/cache logs
- Overall CI duration improvements

## Troubleshooting

If cache isn't being reused:
1. Verify all jobs use same `cache-name`
2. Check that `needs: precompile` is specified
3. Ensure Julia versions match between precompile and test jobs
4. Look for cache key conflicts in job matrix