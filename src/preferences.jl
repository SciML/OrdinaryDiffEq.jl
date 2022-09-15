"""
    set_snoop_precompile!(include_snoop_precompile::Bool = true, export_prefs::Bool = false)

Enable or disable the snoop-precompilation by passing the appropriate boolean value.
By default it is enabled aka `true`.
These preferences are stored in a LocalPrefernces.toml.
Instead, to store it in (OrdinaryDiffEq)Project.toml, set `export_prefs = true`
"""
function set_snoop_precompile!(include_snoop_precompile::Bool = true;
                               export_prefs::Bool = false)
    set_preferences!(@__MODULE__,
                     "include_snoop_precompile" => include_snoop_precompile;
                     export_prefs,
                     force = true)
end

const INCLUDE_SNOOP_PRECOMPILE = load_preference(@__MODULE__,
                                                 "include_snoop_precompile",
                                                 true)
