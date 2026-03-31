module DiffEqBaseGeneralizedGeneratedExt

using DiffEqBase
using GeneralizedGenerated

function SciMLBase.numargs(::GeneralizedGenerated.RuntimeFn{Args}) where {Args}
    return GeneralizedGenerated.from_type(Args) |> length
end

end
