# Repository guidance

- DiffEqBase and OrdinaryDiffEqCore have independent releases. When Core relies on a
  DiffEqBase API or return contract, version that API appropriately and require the
  matching DiffEqBase version in Core's compat lower bound.
- Monorepo `[sources]` entries select local packages. Validate compatibility with
  minimum supported dependency versions, or exclude incompatible versions through
  compat bounds; tests using the checkout do not cover those combinations.
