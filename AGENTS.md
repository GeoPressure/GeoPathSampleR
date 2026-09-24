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

- Performance is a priority.
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

## Error handling

- Do not introduce `tryCatch()` or defensive error handling.
- Let errors occur naturally unless they are handled at a public entry point.

## Silent failure risks

Intervene only when there is a clear risk of unintended recycling, dimension dropping, implicit coercion, or invalid numerical operations. Otherwise, assume correct usage.

## Rewriting vs patching

- Default to a minimal patch.
- If code is clearly inconsistent, inefficient, or unclear, rewriting the function is allowed, but explicitly warn the user.

## Output invariants

Edits must preserve classes, dimensions, ordering, column names, and attributes. Do not silently change them. Explain any necessary trade-off before making it.

## Data structures

- Prefer matrix and array operations when performance matters.
- Avoid unnecessary conversion between matrices and data frames.
- Preserve existing data structures.

## Tests and documentation

- Do not add or modify tests, roxygen documentation, README files, or vignettes unless requested.
- Match validation to the change: source changed functions, run focused tests, and avoid long analyses unless needed.

## Generated files

- Edit roxygen comments in `R/`, never `man/*.Rd` or `NAMESPACE`; run `roxygen2::roxygenise()` and commit the regenerated files.
- Edit `README.Rmd`, never `README.md`; run `devtools::build_readme()` and commit both.
- After editing the vignette or documentation, rebuild with `pkgdown::build_site()` or render the vignette to confirm it runs.
- Keep the software title, version, and DOI consistent across `DESCRIPTION`, `CITATION.cff`, `inst/CITATION`, `codemeta.json`, and `README.Rmd`.

## Mandatory CI checks

Before committing or pushing agent-made changes, run from the repository root:

```sh
jarl check .
air format . --check
```

Both commands must pass without warnings or errors introduced by the changes. If formatting fails, run `air format .`, review the resulting changes, and rerun both checks. Never defer failures to CI.

## Development, PR, and release workflow

### Branch flow

- Develop on `dev`: make focused commits, run the mandatory checks, and push to `origin/dev`.
- Open a pull request from `dev` to `main`. Merge `main` into `dev` and resolve conflicts before the PR is merged.
- Use a draft PR while release metadata or release checks are incomplete. Mark it ready only when the final package version, NEWS entry, and required GitHub checks are complete.

### One canonical release block

- For every release, write one Markdown release block first. It is the source of truth and must be copied **verbatim** to:
  1. the new top section of `NEWS.md`;
  2. the pull-request body; and
  3. the GitHub Release body created for the version tag.
- Do not shorten, paraphrase, reorder, or add items independently in any of those three places.
- Use the release version as the PR title (for example, `v0.2.0`) and as the top NEWS heading.
- Keep the heading, subsection headings, bullet text, Markdown links, and the full-changelog comparison link identical in all three copies.

```md
# GeoPathSampleR vX.Y.Z

## Main

- [Describe the principal user-facing change](https://github.com/GeoPressure/GeoPathSampleR/commit/<sha>).

## Minor

- [Describe a smaller change or fix](https://github.com/GeoPressure/GeoPathSampleR/commit/<sha>).

**Full Changelog**: <https://github.com/GeoPressure/GeoPathSampleR/compare/vX.Y.(Z-1)...vX.Y.Z>
```

### Release checklist

- Set `DESCRIPTION`, `CITATION.cff`, and `codemeta.json` to the final `X.Y.Z` version; do not merge a release with `.9000`.
- Run `cffr::cff_write()` after finalizing `DESCRIPTION`, then review and commit the generated `CITATION.cff` changes. Run `codemetar::write_codemeta()` to update `codemeta.json`.
- Add the canonical release block to `NEWS.md` before opening the PR, then paste that exact block into the PR body.
- Resolve all `R CMD check` warnings and release-relevant notes, and confirm the PR's GitHub Actions matrix is green.
- After merging to `main`, create tag `vX.Y.Z` and paste the unchanged canonical release block into the GitHub Release description. The Zenodo integration archives the release and takes its metadata from `CITATION.cff`.

## Package constraints

- Never use `:::` in package code.
- Call internal GeoPathSampleR functions directly; use `GeoPressureR::` only when a namespace qualifier is needed.
- Do not introduce a helper used only once.
- Do not add dependencies, unnecessary validation, unrelated refactors, silent output changes, or unnecessary data conversions.

## Uncertainty

If a requirement is unclear and could change behaviour, ask before proceeding.
