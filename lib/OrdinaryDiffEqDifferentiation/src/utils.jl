# Reimplementation of InteractiveUtils.subtypes
# Used to remove the direct dependency

subtypes(x::Type) = _subtypes_in!(Base.loaded_modules_array(), x)
function _subtypes_in!(mods::Array, x::Type)
    xt = Base.unwrap_unionall(x)
    if !isabstracttype(x) || !isa(xt, DataType)
        # Fast path
        return Type[]
    end
    sts = Vector{Any}()
    while !isempty(mods)
        m = pop!(mods)
        xt = xt::DataType
        for s in names(m, all = true)
            if isdefined(m, s) && !Base.isdeprecated(m, s)
                t = getfield(m, s)
                dt = isa(t, UnionAll) ? Base.unwrap_unionall(t) : t
                if isa(dt, DataType)
                    if dt.name.name === s && dt.name.module == m && supertype(dt).name == xt.name
                        ti = typeintersect(t, x)
                        ti != Base.Bottom && push!(sts, ti)
                    end
                elseif isa(t, Module) && nameof(t) === s && parentmodule(t) === m && t !== m
                    t === Base || push!(mods, t) # exclude Base, since it also parented by Main
                end
            end
        end
    end
    return permute!(sts, sortperm(map(string, sts)))
end