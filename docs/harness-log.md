# Harness log

Dated notes on what tripped us up and what changed in the harness (make check, hooks, agents,
workflow). Newest last.

## 2026-11-05: Harness added

- Added `make check` (ruff, pytest, smoke run, em-dash check, `bash -n`), AGENTS.md, and the
  `docs/` records. Design decisions moved out of CLAUDE.md into `docs/decisions.md`.
- Tripped: an inline `rm -rf $tmpdir` in the Makefile was blocked by the Claude Code safety
  check. The smoke run now lives in `tests/run_smoke.py`, which uses a self-cleaning temp dir.
- Tripped: writing the Python escape for an em dash through the editor produced the real
  character in a test file, which the new style gate caught. Code that needs the character now uses `chr(0x2014)`.
- Tripped: `make check` takes about 65 s, not 30 s; building the smoke dataset dominates.
- Blocked: pushing to GitHub fails since 2026-11-03 because the Claude GitHub App has no access
  to `ShivPanndey/glm-research`. Commits are local until access is fixed.
