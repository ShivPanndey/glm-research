# `make check` is the one gate to run before every commit and push.
PY ?= python
SHELL_SCRIPTS = $(shell git ls-files '*.sh' '*.sbatch')

.PHONY: check lint test smoke style shellcheck

check: lint test smoke style shellcheck
	@echo "make check: all gates passed"

lint:
	ruff check .

test:
	$(PY) -m pytest -q tests

smoke:
	$(PY) tests/run_smoke.py

style:
	$(PY) scripts/check_style.py

shellcheck:
	@for f in $(SHELL_SCRIPTS); do bash -n "$$f" || exit 1; done; echo "shell syntax ok"
