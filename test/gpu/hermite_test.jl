using ComponentArrays, CUDA, Adapt, RecursiveArrayTools, FastBroadcast, FillArrays,
      OrdinaryDiffEq, Test

a = ComponentArray((a = rand(Float32, 5, 5), b = rand(Float32, 5, 5)))
a = adapt(CuArray, a)
pa = ArrayPartition(a)
pb = deepcopy(pa)
pc = deepcopy(pa)
pd = deepcopy(pa)
pe = deepcopy(pa)
k = [pd, pe]
t = FillArrays.Trues(length(pa))

OrdinaryDiffEq.hermite_interpolant!(pa, 0.1, 0.2, pb, pc, k, nothing, Val{0}, t) # if this doesn't error we're good

@test pa.x[1] != pb.x[1]
