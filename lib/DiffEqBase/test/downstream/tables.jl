using OrdinaryDiffEq, DataFrames, Test, SymbolicIndexingInterface
f_2dlinear = (du, u, p, t) -> du .= 1.01u;
prob = ODEProblem(f_2dlinear, rand(2, 2), (0.0, 1.0));
sol1 = solve(prob, Euler(); dt = 1 // 2^(4));
df = DataFrame(sol1)
@test names(df) == ["timestamp", "value1", "value2", "value3", "value4"]

prob = ODEProblem(
    ODEFunction(
        f_2dlinear, sys = SymbolicIndexingInterface.SymbolCache([:a, :b, :c, :d], [], :t)
    ),
    rand(2, 2),
    (0.0, 1.0)
);
sol2 = solve(prob, Euler(); dt = 1 // 2^(4));
df = DataFrame(sol2)
@test names(df) == ["timestamp", "a", "b", "c", "d"]
