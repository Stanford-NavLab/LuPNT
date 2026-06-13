# GitHub Actions

This directory contains the GitHub Actions workflows for LuPNT. The active CI
workflows are built around Pixi so the same named tasks can be run locally and
in CI.

Most CI workflows run on pushes and pull requests targeting the `development`
branch. Release publishing is tag-based, and the Windows workflow is currently a
manual status workflow.

## Workflow Summary

| Workflow | File | Trigger | What it does |
| --- | --- | --- | --- |
| Style | `style.yml` | Push and pull request to `development` | Installs the Pixi environment and runs `pixi run format-check` against the pull request base or previous commit. |
| Ubuntu | `ubuntu.yml` | Push and pull request to `development` | Restores cached data and CPM sources, then runs the C++ CI suite with `pixi run test-cpp-ci` on Ubuntu. |
| MacOS | `macos.yml` | Push and pull request to `development` | Runs the same C++ CI suite as Ubuntu with `pixi run test-cpp-ci` on macOS. |
| Python | `python.yml` | Push and pull request to `development` | Runs `pixi run test-py` on both Ubuntu and macOS. |
| Examples | `examples.yml` | Push and pull request to `development` | Builds the example programs with `pixi run build-examples`. |
| Install | `install.yml` | Push and pull request to `development` | Runs the install smoke test through `pixi run install-test`. |
| Docs | `docs.yaml` | Push and pull request to `development` | Builds the documentation with `pixi run build-docs`; on pushes, uploads the generated site to GitHub Pages. |
| PyPI | `pypi.yml` | Version tags matching `v*` and manual dispatch | Builds a source distribution and publishes it to PyPI on version tag pushes. |
| Windows | `windows.yml` | Manual dispatch | Reports that Windows CI is paused until `win-64` Pixi support is enabled and validated. |

## Shared CI Behavior

- Workflows use `prefix-dev/setup-pixi` with the lockfile in frozen mode.
- C++ workflows cache CPM sources through `CPM_SOURCE_CACHE`.
- Test workflows cache `data/LuPNT_data` and verify key files exist before
  reusing the cache, which avoids silently using partial data downloads.
- Pixi activation provides environment variables such as `LUPNT_DATA_PATH`.
- CI task definitions live in `pixi.toml`; update this README when task names or
  workflow responsibilities change.

## Documentation Deployment

The docs workflow builds the Sphinx documentation into `build/docs` and deploys
that directory with the official GitHub Pages artifact flow:

1. `actions/configure-pages`
2. `actions/upload-pages-artifact`
3. `actions/deploy-pages`

This avoids pushing directly to the protected `gh-pages` branch. The docs build
also removes Sphinx doctree cache files from the published artifact, and limits
the docs build parallelism with `LUPNT_DOCS_JOBS=4`.

## Running Checks Locally

Useful local equivalents:

```bash
pixi install --frozen
pixi run format-check
pixi run test-cpp-ci
pixi run test-py
pixi run build-examples
pixi run install-test
pixi run build-docs
```

For release validation, use a version tag such as `v0.1.0`. PyPI publishing
requires the repository to be configured for PyPI trusted publishing.
