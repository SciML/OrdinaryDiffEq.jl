using OrdinaryDiffEqExtrapolation
using SciMLBase: ODEFunction
using Test

const ODEX = OrdinaryDiffEqExtrapolation

@testset "Local dependency-independent utilities" begin
    @test ODEX._unwrap_val(Val(true))
    @test !ODEX._unwrap_val(Val(false))
    @test ODEX._unwrap_val(nothing) === nothing

    vector = [1.0, 2.0]
    matrix = [1.0 2.0; 3.0 4.0]
    @test ODEX._vec(1.0) === 1.0
    @test ODEX._vec(vector) === vector
    @test ODEX._vec(matrix) == vec(matrix)
    @test ODEX._reshape(1.0, (1,)) === 1.0
    @test ODEX._reshape(vector, (1, 2)) == reshape(vector, 1, 2)

    f = ODEFunction((u, p, t) -> u)
    f_with_Wfact = ODEFunction((u, p, t) -> u; Wfact = (u, p, dtgamma, t) -> u)
    @test !ODEX._has_Wfact(f)
    @test ODEX._has_Wfact(f_with_Wfact)

    @test ODEX._controller_scalar_type((;)) === Float64
    @test ODEX._controller_scalar_type((; qmax = 2)) === Float64
    @test ODEX._controller_scalar_type((; qmax = Float32(2))) === Float32

    @test ODEX._thread_storage_size() ==
        Threads.nthreads(:default) + Threads.nthreads(:interactive)
end

@testset "Nested threaded fallback keeps caller migration state" begin
    for throws in (false, true)
        before = fill(true, 2)
        after = fill(false, 2)
        stable = fill(false, 2)
        caught = fill(false, 2)
        Threads.@threads :dynamic for i in 1:2
            task = current_task()
            before[i] = task.sticky
            ids = Int[]
            try
                ODEX._threaded_foreach(true, 1:4) do j
                    push!(ids, Threads.threadid())
                    yield()
                    push!(ids, Threads.threadid())
                    throws && j == 3 && error("body failed")
                end
            catch
                caught[i] = true
            end
            after[i] = task.sticky
            stable[i] = allequal(ids)
        end
        @test before == after == [false, false]
        @test all(stable)
        @test caught == [throws, throws]
    end
end
