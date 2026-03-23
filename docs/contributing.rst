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

   $ make lint        # formatting, import order, static analysis (pre-commit)
   $ make test        # pytest with coverage data collection
   $ make coverage    # combine coverage files and print report
   $ make typecheck   # runtime type checking via typeguard
   $ make doctest     # validate inline docstring examples

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

.. _sphinx: https://www.sphinx-doc.org/en/master/

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
