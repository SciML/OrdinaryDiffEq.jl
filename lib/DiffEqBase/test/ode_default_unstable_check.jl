using Test, RecursiveArrayTools, StaticArrays, SparseArrays

using DiffEqBase: NAN_CHECK

@test !NAN_CHECK(3.0 + 4.0im)
@test NAN_CHECK(NaN)

u1 = ones(3)
@test !NAN_CHECK(u1)
u1′ = copy(u1)
u1′[2] = NaN
@test NAN_CHECK(u1′)

u2 = [SA[1.0 1.0; 1.0 1.0] for i in 1:3]
@test !NAN_CHECK(u2)
u2′ = copy(u2)
u2′[2] = SA[1.0 NaN; 1.0 1.0]
@test NAN_CHECK(u2′)

u3 = VectorOfArray([ones(5), ones(5)])
@test !NAN_CHECK(u3)
u3′ = recursivecopy(u3)
u3′[3, 2] = NaN
@test NAN_CHECK(u3′)

u4 = ArrayPartition(u1, u2, u3)
@test !NAN_CHECK(u4)
u4_1 = ArrayPartition(u1′, u2, u3)
@test NAN_CHECK(u4_1)
u4_2 = ArrayPartition(u1, u2′, u3)
@test NAN_CHECK(u4_2)
u4_3 = ArrayPartition(u1, u2, u3′)
@test NAN_CHECK(u4_3)

@test !NAN_CHECK(ArrayPartition(u4, u4))
@test NAN_CHECK(ArrayPartition(u4, u4_1))
@test NAN_CHECK(ArrayPartition(u4, u4_2))
@test NAN_CHECK(ArrayPartition(u4, u4_3))

u5 = spzeros(1, 1)
@test !NAN_CHECK(u5)
u5[1, 1] = NaN
@test NAN_CHECK(u5)
