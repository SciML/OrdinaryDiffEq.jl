"""
    autoplot(wp_set; families=nothing, reference_tags=nothing, best_n=2)

Generate a collection of WorkPrecisionSets for standard comparison plots.

Returns a `Dict{String, WorkPrecisionSet}` with:
- `"family_<tag>"` for each family
- `"best_of_families"` with top methods per family
- `"all"` containing the original set
"""
function autoplot(wp_set::WorkPrecisionSet;
        families = nothing,
        reference_tags = nothing,
        best_n::Int = 2)
    results = Dict{String, WorkPrecisionSet}()

    # Determine families
    if families === nothing
        families = _auto_detect_families(wp_set)
    end

    # Per-family plots
    for fam in families
        subset = filter_by_tags(wp_set, fam)
        if reference_tags !== nothing
            ref_tags = reference_tags isa Symbol ? (reference_tags,) : Tuple(reference_tags)
            ref = filter_by_tags(wp_set, ref_tags...)
            # Deduplicate
            seen = Set(wp.name for wp in subset.wps)
            extra_wps = [wp for wp in ref.wps if wp.name ∉ seen]
            if !isempty(extra_wps)
                subset = merge_wp_sets(subset,
                    WorkPrecisionSet(extra_wps, length(extra_wps),
                        wp_set.abstols, wp_set.reltols, wp_set.prob,
                        wp_set.setups, [wp.name for wp in extra_wps],
                        wp_set.error_estimate, wp_set.numruns))
            end
        end
        length(subset) > 0 && (results["family_$(fam)"] = subset)
    end

    # Best-of-families
    if !isempty(families)
        best_wps = WorkPrecision[]
        seen = Set{String}()
        for fam in families
            subset = filter_by_tags(wp_set, fam)
            length(subset) == 0 && continue
            ranked = sort(collect(1:length(subset.wps)),
                by = i -> let t = filter(!isnan, subset.wps[i].times)
                    isempty(t) ? Inf : mean(t)
                end)
            for idx in ranked[1:min(best_n, length(ranked))]
                wp = subset.wps[idx]
                if wp.name ∉ seen
                    push!(seen, wp.name)
                    push!(best_wps, wp)
                end
            end
        end
        if !isempty(best_wps)
            names = [wp.name for wp in best_wps]
            results["best_of_families"] = WorkPrecisionSet(
                best_wps, length(best_wps), wp_set.abstols, wp_set.reltols,
                wp_set.prob, wp_set.setups, names, wp_set.error_estimate, wp_set.numruns)
        end
    end

    results["all"] = wp_set
    return results
end

function _auto_detect_families(wp_set::WorkPrecisionSet)
    all_tags = unique_tags(wp_set)
    n_methods = length(wp_set.wps)
    n_methods == 0 && return Symbol[]
    # Filter out tags that appear on >80% of methods (likely category tags)
    return filter(all_tags) do tag
        count(wp -> tag in wp.tags, wp_set.wps) <= 0.8 * n_methods
    end
end
