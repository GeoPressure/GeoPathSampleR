# AGENTS.md

## Scope

- Applies to the whole repository.
- Covers code edits, tests, and documentation.
- Default to minimal, surgical modifications unless explicitly asked otherwise.

## Language and dependencies

- Use R (>= 4.1) and the native pipe `|>`.
- Do not introduce new dependencies.
- Use existing dependencies where appropriate.
- Use `glue::glue()` for string interpolation and avoid `paste()` unless performance requires it.

## General principles

- Assume inputs are valid and well-formed.
- Do not add defensive programming or guards unless explicitly requested.
- Prefer the shortest correct implementation.
- Avoid unnecessary line breaks, verbosity, and variables.
- Introduce a variable only when it prevents recomputation or is required for correctness or performance.

## Performance

- Prefer vectorised operations, matrix or array operations, `rowSums`, `colSums`, `sweep`, `%*%`, and `outer` where they fit.
- Prefer efficient primitives over `apply()` when possible. Use `apply()` or `Map()` only when they remain clear and appropriate.
- Avoid unnecessary copies, repeated computation, and hidden coercions.
- Do not trade memory for speed, or speed for memory, unless explicitly requested.

## Code style and structure

### Pipes

- Prefer compact pipe chains: `x |> f() |> g() |> h()`.
- Do not introduce intermediate variables unless they improve correctness, performance, or inspectability.

### Comments

- Add minimal comments per logical section.
- Describe intent rather than obvious syntax.
- Use short section headers where helpful.

### Function edits

- Make the smallest possible patch.
- Do not refactor or reorder unrelated code.
- Keep reusable functions in `R/`; keep executable analysis code in scripts or vignettes.

## Validation and numerical behaviour

- Validate user-facing inputs at public entry points or dedicated validation
  helpers invoked immediately from them.
- Keep computational helpers free of speculative validation and error handling.
- Preserve numerical behaviour unless there is a clear risk such as `log(0)`, division by zero, or unstable normalisation.
- When adding a numerical safeguard, keep it minimal and explain why it is needed.

## Output invariants

Edits must preserve classes, dimensions, ordering, column names, and attributes. Do not silently change them. Explain any necessary trade-off before making it.

## Data structures

- Prefer matrix and array operations when performance matters.
- Avoid unnecessary conversion between matrices and data frames.
- Preserve existing data structures.

## Tests and documentation

- Do not add or modify tests, roxygen documentation, README files, or vignettes unless requested.
- Match validation to the change: source changed functions, run focused tests, and avoid long analyses unless needed.

## Mandatory CI checks

Before committing or pushing agent-made changes, run from the repository root:

```sh
jarl check .
air format . --check
```

Both commands must pass without warnings or errors introduced by the changes. If formatting fails, run `air format .`, review the resulting changes, and rerun both checks. Never defer failures to CI.

## Package constraints

- Never use `:::` in package code.
- Call internal GeoPathSampleR functions directly; use `GeoPressureR::` only when a namespace qualifier is needed.
- Do not introduce a helper used only once.
- Do not add dependencies, unnecessary validation, unrelated refactors, silent output changes, or unnecessary data conversions.

## Uncertainty

If a requirement is unclear and could change behaviour, ask before proceeding.
