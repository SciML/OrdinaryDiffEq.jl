# OrdinaryDiffEq.jl QA Testing Setup

This document describes the comprehensive QA testing infrastructure that has been set up for OrdinaryDiffEq.jl to verify non-allocating behavior and type stability of solvers during stepping operations.

## Overview

The QA testing setup includes:

1. **Allocation Tests**: Verify that `step!` operations don't allocate memory during stepping
2. **Type Stability Tests**: Use JET.jl to verify type stability of solver operations
3. **Systematic Testing**: Test all solver sublibraries with appropriate test problems
4. **@test_broken Framework**: Track solvers that are currently allocating or type-unstable

## Key Findings

### Currently Non-Allocating Solvers (during basic testing)
Based on initial testing, these solvers showed allocation-free behavior in simple scenarios:
- **Explicit RK**: RK4, BS3, DP5, Vern6, Vern7, Vern8, Vern9
- **SSPRK**: SSPRK43
- **Low-order**: Euler (with fixed timestep)

### Currently Allocating Solvers
These solvers are currently allocating during `step!` operations and marked with `@test_broken`:
- **Tsit5**: 16 bytes per step (unexpected - needs investigation)
- **BDF Methods**: QNDF, ABDF2 (expected for implicit methods)
- **Rosenbrock**: Rodas4, Rodas5 (expected for implicit methods)  
- **SSPRK**: SSPRK22, SSPRK33
- **SDIRK**: ImplicitEuler, KenCarp4, SDIRK2, TRBDF2

### Type Stability Status
- **All tested solvers** currently show type instabilities and are marked with `@test_broken`
- Type instabilities are detected in initialization, stepping, and full solve operations

## Test File Structure

The following test files have been created in the respective sublibrary directories:

```
lib/
├── OrdinaryDiffEqTsit5/test/
│   ├── allocation_tests.jl    # Allocation tests for Tsit5
│   └── jet_tests.jl          # JET type stability tests
├── OrdinaryDiffEqExplicitRK/test/
│   └── allocation_tests.jl    # Tests for RK4, BS3, DP5
├── OrdinaryDiffEqHighOrderRK/test/
│   └── allocation_tests.jl    # Tests for Vern6-9
├── OrdinaryDiffEqLowOrderRK/test/
│   └── allocation_tests.jl    # Tests for Euler, RK4, etc.
├── OrdinaryDiffEqSSPRK/test/
│   └── allocation_tests.jl    # Tests for SSPRK methods
├── OrdinaryDiffEqBDF/test/
│   └── allocation_tests.jl    # Tests for BDF methods (@test_broken)
└── OrdinaryDiffEqRosenbrock/test/
    └── allocation_tests.jl    # Tests for Rosenbrock methods (@test_broken)
```

## Test Framework Components

### 1. Allocation Testing Pattern

```julia
@testset "Solver Allocation Tests" begin
    # Test problem setup
    function simple_system!(du, u, p, t)
        du[1] = -0.5 * u[1]
        du[2] = -1.5 * u[2]
    end
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    # Test allocation-free solvers
    allocation_free_solvers = [RK4(), BS3(), DP5()]
    
    for solver in allocation_free_solvers
        @testset "$(typeof(solver)) allocation test" begin
            integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
            step!(integrator)  # Setup step may allocate
            
            # Test subsequent steps for zero allocations
            for i in 2:10
                if integrator.t >= integrator.sol.prob.tspan[2]
                    break
                end
                alloc = @allocated step!(integrator)
                @test alloc == 0  # Should be allocation-free
            end
        end
    end
    
    # Test currently allocating solvers with @test_broken
    allocating_solvers = [Tsit5(), ImplicitEuler()]
    
    for solver in allocating_solvers
        @testset "$(typeof(solver)) allocation test (broken)" begin
            integrator = init(prob, solver, save_everystep=false, abstol=1e-6, reltol=1e-6)
            step!(integrator)  # Setup step may allocate
            
            for i in 2:6
                if integrator.t >= integrator.sol.prob.tspan[2]
                    break
                end
                alloc = @allocated step!(integrator)
                @test_broken alloc == 0  # Currently allocating, should be fixed
            end
        end
    end
end
```

### 2. JET Type Stability Testing Pattern

```julia
using JET

@testset "JET Type Stability Tests" begin
    prob = ODEProblem(simple_system!, [1.0, 1.0], (0.0, 1.0))
    
    @testset "Solver Initialization Type Stability" begin
        @test_opt broken=true init(prob, Tsit5(), save_everystep=false, abstol=1e-6, reltol=1e-6)
    end
    
    @testset "Step Operation Type Stability" begin
        integrator = init(prob, Tsit5(), save_everystep=false, abstol=1e-6, reltol=1e-6)
        @test_opt broken=true step!(integrator)
    end
    
    @testset "Full Solve Type Stability" begin
        @test_opt broken=true solve(prob, Tsit5(), save_everystep=false, abstol=1e-6, reltol=1e-6)
    end
end
```

## Running the Tests

### Individual Sublibrary Tests
```bash
# Test Tsit5 allocations
cd lib/OrdinaryDiffEqTsit5
julia --project=../.. -e 'include("test/allocation_tests.jl")'

# Test Tsit5 type stability  
julia --project=../.. -e 'include("test/jet_tests.jl")'
```

### Comprehensive Testing
Run the comprehensive test script to analyze all solvers:
```bash
julia --project test_comprehensive_qa.jl
```

## Key Dependencies Added

The following testing dependencies have been added to Project.toml:
- `AllocCheck.jl`: Static analysis for allocation-free code verification
- `JET.jl`: Static analysis for type stability verification

## Usage for Development

### For Solver Developers

1. **When implementing new solvers**: Add allocation tests following the patterns above
2. **When fixing allocation issues**: Move solvers from `@test_broken` to regular `@test` assertions
3. **When fixing type stability**: Remove `broken=true` from JET tests

### Adding New Tests

1. Create `allocation_tests.jl` in the sublibrary's `test/` directory
2. Follow the established patterns for test structure
3. Use `@test_broken` for currently failing tests
4. Include both allocation and type stability tests

## Expected Development Workflow

1. **Current State**: Most solvers have allocation or type stability issues (marked `@test_broken`)
2. **Development Goal**: Fix underlying allocation and type stability issues
3. **Progress Tracking**: Convert `@test_broken` to `@test` as issues are resolved
4. **Verification**: Use these tests to ensure solvers remain allocation-free after fixes

## Test Problem Selection

- **Explicit solvers**: Use vector ODEs like Lorenz or simple linear systems
- **Implicit/stiff solvers**: Use scalar linear problems or stiff ODEs
- **Fixed timestep methods**: Specify `dt` parameter and set `adaptive=false`

## Integration with CI

These tests can be integrated into the existing CI pipeline:

1. Add to sublibrary `runtests.jl` files
2. Run during PR testing to catch regressions
3. Track progress on allocation-free and type-stable solver development

## Future Extensions

1. **Performance benchmarking**: Extend to measure solve times alongside allocations
2. **Memory profiling**: Add detailed memory usage analysis
3. **Regression testing**: Ensure fixes don't break other functionality
4. **Documentation**: Generate reports showing solver performance characteristics

This QA testing infrastructure provides a solid foundation for systematically improving the performance characteristics of OrdinaryDiffEq.jl solvers.