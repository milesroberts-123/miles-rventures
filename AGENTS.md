# AGENTS.md

## Project

An R package of convenience functions for use across RStudio projects (per `README.md`).

- Package name: `rventures` (see `DESCRIPTION`).
- Roxygen2 is used; after editing `R/` files, run `devtools::document()` to regenerate `NAMESPACE` and `man/`.
- `R/` files are organized by module, not one file per function: `plotting.R`, `file_management.R`, `statistics.R`, `popgen.R`. Add new functions to the matching module.

## Toolchain

- Language: R (RStudio workflow implied by `.gitignore`).
- License: MIT, copyright Miles Roberts.

## Git

- Branch: `main`; remote `origin` → `github.com/milesroberts-123/miles-rventures.git`.
- Never push to `main` without asking first.
