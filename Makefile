OKTOBERFEST_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
# include $(OKTOBERFEST_DIR)MakefileShared # sourcing this is not needed because it is already sourced in MakefileCe
include $(OKTOBERFEST_DIR)MakefileShared

dependencies:
	git config --global credential.helper cache

registry:
	docker login gitlab.lrz.de:5005
	docker build -t gitlab.lrz.de:5005/proteomics/github/oktoberfest .
	docker push gitlab.lrz.de:5005/proteomics/github/oktoberfest

jump:
	$(DOCKER_CMD) \
		$(IMAGE) bash

# --no-cache
build: dependencies
	git describe --long --dirty --always > hash.file
	docker build -f Dockerfile -t $(IMAGE) . || (exit 1)

bootstrap: DATA=/root/data
bootstrap:
	bash -c "cp /root/Makefile* $(LOCAL_DIR)"

run_oktoberfest: rm_err_file
	$(DOCKER_CMD) \
		$(IMAGE) python3 -u -m oktoberfest --config_path $(LOCAL_DIR)/config.json || (echo "2" > $(DATA)err.out; exit 2)

compress: run_oktoberfest
	zip -j -r -9 "$(DATA)/results.zip" "$(DATA)/results/" || (echo "3" > $(DATA)err.out; exit 3)

all: compress


run_local:
	python3 -u -m oktoberfest --config_path $(DATA)/config.json

clean_data_folder:
	bash -c "rm -rf $(DATA)/{proc,msms,results,mzML,msms.prosit,err.out,results.zip}"


# ── Developer / CI targets ──────────────────────────────────────────────────
# These mirror exactly what CI runs. Use `make install` once to set up, then
# any target below can be run locally or called from a workflow step.
#
# Usage:
#   make install     – install project + all dev deps into Poetry virtualenv
#   make check       – run all quality checks in CI order (lint → test → coverage → typecheck → doctest)
#   make lint        – run all pre-commit hooks (format, lint, import order)
#   make test        – run tests and collect coverage data
#   make coverage    – combine .coverage.* files and print report
#   make typecheck   – runtime type checking via typeguard
#   make doctest     – validate inline docstring examples
#   make docs        – build HTML documentation with Sphinx
#   make dist        – build source and wheel distributions
# ────────────────────────────────────────────────────────────────────────────
.PHONY: install lint test coverage typecheck doctest docs dist check

check: lint test coverage typecheck doctest ## Run all quality checks in CI order

install: ## Install project and all dev dependencies
	poetry install

lint: ## Run all pre-commit hooks (formatting, linting, import order)
	pre-commit run --all-files

test: ## Run test suite and collect coverage data
	poetry run coverage run --parallel -m pytest tests/

coverage: ## Combine coverage files and print report
	poetry run coverage combine
	poetry run coverage report -i

typecheck: ## Runtime type checking with typeguard
	poetry run pytest --typeguard-packages=oktoberfest tests/

doctest: ## Validate inline docstring examples with xdoctest
	poetry run python -m xdoctest oktoberfest all

docs: ## Build HTML documentation with Sphinx
	sphinx-build -b html docs docs/_build

dist: ## Build source and wheel distributions
	poetry build --ansi
