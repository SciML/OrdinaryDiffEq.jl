"""
    autoplot(wp_set::WorkPrecisionSet; families = nothing, reference_tags = nothing,
        best_n = 2) -> Dict{String, WorkPrecisionSet}

Split one tagged [`WorkPrecisionSet`](@ref) into the subsets a benchmark page usually
plots, without re-solving anything:

  - `"family_<tag>"` for each family in `families`, holding that family's entries plus
    any entry matching `reference_tags`
  - `"best_of_families"`, the [`best_of_families`](@ref) selection of the `best_n` best
    entries per family, again alongside the reference entries
  - `"all"`, the full input set

Each value plots directly with `plot(subset)`. When `families` is not given it is taken
from the tags in use, dropping tags carried by more than 80% of the entries, since those
describe the benchmark as a whole (`:stiff`, `:nonstiff`) rather than a family.
"""
function autoplot(
        wp_set::WorkPrecisionSet;
        families = nothing,
        reference_tags = nothing,
        best_n::Int = 2
    )
    families === nothing && (families = _auto_detect_families(wp_set))
    ref_indices = reference_tags === nothing ? Int[] :
        findall(wp -> _has_any_tag(wp, _as_tags(reference_tags)), wp_set.wps)

    results = Dict{String, WorkPrecisionSet}()
    present = Symbol[]
    for family in families
        indices = findall(wp -> family in wp.tags, wp_set.wps)
        isempty(indices) && continue
        push!(present, family)
        results["family_$(family)"] = _subset_wps(wp_set, sort!(union(indices, ref_indices)))
    end

    if !isempty(present)
        best = _best_family_indices(wp_set, present, best_n, :area)
        results["best_of_families"] = _subset_wps(wp_set, sort!(union(best, ref_indices)))
    end

    results["all"] = wp_set
    return results
end

function _auto_detect_families(wp_set::WorkPrecisionSet)
    n_methods = length(wp_set.wps)
    n_methods == 0 && return Symbol[]
    return filter(unique_tags(wp_set)) do tag
        count(wp -> tag in wp.tags, wp_set.wps) <= 0.8 * n_methods
    end
end
