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

@testset "Threaded foreach gives each chunk its own slot" begin
    for range in (1:2, 0:1, 0:2, 0:5), nested in (false, true)
        nchunks = min(length(range), Threads.nthreads())
        outer = nested ? 2 : 1
        slots = [zeros(Int, length(range)) for _ in 1:outer]
        visits = [zeros(Int, length(range)) for _ in 1:outer]
        busy = [fill(false, nchunks) for _ in 1:outer]
        clash = fill(false, outer)
        Threads.@threads :dynamic for k in 1:outer
            ODEX._threaded_foreach(true, range) do i, slot
                clash[k] |= busy[k][slot]
                busy[k][slot] = true
                yield()
                slots[k][i - first(range) + 1] = slot
                visits[k][i - first(range) + 1] += 1
                busy[k][slot] = false
            end
        end
        for k in 1:outer
            sizes = [count(==(c), slots[k]) for c in 1:nchunks]
            @test all(==(1), visits[k])
            @test issorted(slots[k])
            @test all(>(0), sizes) && maximum(sizes) - minimum(sizes) <= 1
            @test !clash[k]
        end
    end
end
