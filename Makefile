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
#   make lint        – run pre-commit hooks (formatting, linting, security checks)
#   make format      – format code with ruff
#   make test        – run tests and collect coverage data
#   make coverage    – print coverage report and export as XML for CI upload
#   make typecheck   – runtime type checking via typeguard
#   make doctest     – validate inline docstring examples
#   make docs        – build HTML documentation with Sphinx
#  make docs-serve  – build docs and serve locally with live reload
#   make check       – run all quality checks in CI order (lint → test → coverage → typecheck → doctest)
#   make dist        – build source and wheel distributions
# ────────────────────────────────────────────────────────────────────────────
.PHONY: install lint format test coverage typecheck doctest docs dist check

check: lint test coverage typecheck doctest ## Run all quality checks in CI order

install: ## Install project with all dev and docs dependencies
	poetry install --with dev --extras docs

lint: ## Run pre-commit hooks (formatting, linting, security checks)
	poetry run pre-commit run --all-files

format: ## Format code with ruff
	poetry run ruff format oktoberfest tests

test: ## Run test suite and collect coverage data
	poetry run coverage run -m pytest tests/unit_tests

coverage: ## Generate coverage report and export as XML
	poetry run coverage report -i
	poetry run coverage xml

typecheck: ## Runtime type checking with typeguard
	poetry run pytest --typeguard-packages=oktoberfest tests/unit_tests

doctest: ## Validate inline docstring examples with xdoctest
	poetry run python -m xdoctest oktoberfest all

docs: ## Build HTML documentation with Sphinx
	poetry run sphinx-build -b html docs docs/_build

docs-serve:  ## build docs and serve locally with live reload
	poetry run sphinx-autobuild docs docs/_build/html --open-browser

dist: ## Build source and wheel distributions
	poetry build --ansi
