# AGENTS.md

## Project

An R package of convenience functions for use across RStudio projects (per `README.md`).

- Package name: `rventures` (see `DESCRIPTION`).
- Roxygen2 is used; after editing `R/` files, run `devtools::document()` to regenerate `NAMESPACE` and `man/`.
- `R/` files are organized by module, not one file per function: `plotting.R`, `file_management.R`, `statistics.R`, `popgen.R`. Add new functions to the matching module.
- Tests use testthat (edition 3). Run all with `devtools::test()`, or one module with `devtools::test(filter = "popgen")`. Test files live in `tests/testthat/`, one per module.

## Toolchain

- Language: R (RStudio workflow implied by `.gitignore`).
- License: MIT, copyright Miles Roberts.

## Git

- Branch: `main`; remote `origin` → `github.com/milesroberts-123/miles-rventures.git`.
- Never push to `main` without asking first.
- Before every push: export the full OpenCode transcript as Markdown (`session-<short-id>.md` in the repo root), compress it with `xz -9` alongside, and commit both files.
