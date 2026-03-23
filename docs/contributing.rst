Contributor Guide
=================

Thank you for your interest in improving this project.
This project is open-source under the `MIT license`_ and
highly welcomes contributions in the form of bug reports, feature requests, and pull requests.

Here is a list of important resources for contributors:

- `Source Code`_
- `Documentation`_
- `Issue Tracker`_
- `Code of Conduct`_

.. _MIT license: https://opensource.org/licenses/MIT
.. _Source Code: https://github.com/wilhelm-lab/oktoberfest
.. _Documentation: https://oktoberfest.readthedocs.io/
.. _Issue Tracker: https://github.com/wilhelm-lab/oktoberfest/issues

How to report a bug
-------------------

Report bugs on the `Issue Tracker`_.


How to request a feature
------------------------

Request features on the `Issue Tracker`_.


How to set up your development environment
------------------------------------------

You need Python 3.10+ and Poetry_:

.. code:: console

    $ pip install poetry

Install the package with all development requirements:

.. code:: console

   $ make install

You can now run an interactive Python session or the command-line interface:

.. code:: console

   $ poetry run python
   $ poetry run oktoberfest

.. _Poetry: https://python-poetry.org/


How to test the project
-----------------------

Run all quality checks in the same order as CI:

.. code:: console

   $ make check

Or run individual steps:

.. code:: console

   $ make lint        # pre-commit hooks (formatting, linting, security checks)
   $ make format      # format code with ruff
   $ make test        # run test suite and collect coverage data
   $ make coverage    # generate coverage report and export as XML
   $ make typecheck   # runtime type checking with typeguard
   $ make doctest     # validate inline docstring examples with xdoctest

Unit tests are located in the ``tests`` directory and are written using the
pytest_ testing framework.

.. _pytest: https://pytest.readthedocs.io/

How to build and view the documentation
---------------------------------------

This project uses Sphinx_ together with several extensions to build the documentation.

Build the HTML documentation from the repository root:

.. code:: console

    $ make docs

The generated static HTML files are in ``docs/_build/html``.
Open ``docs/_build/html/index.html`` in your browser to inspect them.

Build and serve with live reload:

.. code:: console

    $ make docs-serve

This opens your browser and automatically rebuilds and reloads the documentation when files change.

.. _sphinx: https://www.sphinx-doc.org/en/master/

How to make a release
---------------------

Releases are published to PyPI automatically when a GitHub Release is published.
The version string lives only in ``pyproject.toml`` — ``__version__`` is read from
the installed package metadata at runtime.

Release Drafter continuously updates a draft GitHub Release with an accumulated changelog
from merged PR labels and a suggested next version (e.g. ``0.9.1``). It is a changelog
generator — it never modifies any file in the repository.

**Branch model:** ``development`` is the integration branch; ``main`` mirrors exactly
what is published on PyPI. The release tag is created on ``development`` and subsequently
merged into ``main``.

1. **Check the draft release** on GitHub to see the suggested next version (e.g. ``0.10.0``).
   The version is inferred automatically from the labels on merged PRs since the last release.

2. **Bump the version on** ``development``:

   .. code:: console

      $ git checkout development && git pull
      $ poetry version <next-version>   # e.g. poetry version 0.10.0
      $ git add pyproject.toml
      $ git commit -m "bump version to $(poetry version -s)"
      $ git push origin development

3. **Publish the draft release** on GitHub.
   The draft already targets ``development`` (set via ``commitish: development`` in
   ``.github/release-drafter.yml``), so no target branch change is needed.
   Clicking **Publish release** triggers the publish workflow, which:

   - Re-runs the full CI suite as a hard gate.
   - Builds the wheel and sdist with ``poetry build``.
   - Publishes to PyPI via OIDC Trusted Publishing (no secrets required).
   - Creates the tag ``v<next-version>`` on ``development``.

4. **Merge the tagged commit into** ``main`` so that ``main`` reflects the release:

   .. code:: console

      $ git checkout main && git pull
      $ git merge v<next-version> --no-ff -m "release: v<next-version>"
      $ git push origin main

How to submit changes
---------------------

Open a `pull request`_ to submit changes to this project against the ``development`` branch.

Your pull request needs to meet the following guidelines for acceptance:

- All checks in ``make check`` must pass without errors or warnings.
- Include unit tests. This project maintains a high code coverage.
- If your changes add functionality, update the documentation accordingly.

To install pre-commit as a Git hook so checks run automatically on every commit:

.. code:: console

   $ poetry run pre-commit install

It is recommended to open an issue before starting work on anything.
This will allow a chance to talk it over with the owners and validate your approach.

.. _pull request: https://github.com/wilhelm-lab/oktoberfest/pulls
.. _Code of Conduct: CODE_OF_CONDUCT.rst
