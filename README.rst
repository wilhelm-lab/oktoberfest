Oktoberfest: Rescoring, Collision Energy Calibration and Spectral Library Generation for Proteomics
===================================================================================================

|PyPI| |Python Version| |License| |Read the Docs| |CI| |Codecov| |pre-commit| |Ruff|

.. |PyPI| image:: https://img.shields.io/pypi/v/oktoberfest.svg
   :target: https://pypi.org/project/oktoberfest/
   :alt: PyPI
.. |Python Version| image:: https://img.shields.io/pypi/pyversions/oktoberfest
   :target: https://pypi.org/project/oktoberfest
   :alt: Python Version
.. |License| image:: https://img.shields.io/github/license/wilhelm-lab/oktoberfest
   :target: https://opensource.org/licenses/MIT
   :alt: License
.. |Read the Docs| image:: https://img.shields.io/readthedocs/oktoberfest/latest.svg?label=Read%20the%20Docs
   :target: https://oktoberfest.readthedocs.io/
   :alt: Read the documentation at https://oktoberfest.readthedocs.io/
.. |CI| image:: https://github.com/wilhelm-lab/oktoberfest/workflows/CI/badge.svg
   :target: https://github.com/wilhelm-lab/oktoberfest/actions?workflow=CI
   :alt: CI Status
.. |Codecov| image:: https://codecov.io/gh/wilhelm-lab/oktoberfest/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/wilhelm-lab/oktoberfest
   :alt: Codecov
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit
.. |Ruff| image:: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
   :target: https://github.com/astral-sh/ruff
   :alt: Ruff

Oktoberfest is a python tool for collision energy calibration, rescoring search results and generating spectral libraries for proteomics research within the Prosit ecosystem. It offers an end to end pipeline that takes search results, predicts peptide properties using koina, plots summaries and quality control figures and performs FDR estimation with either mokapot or percolator.

Documentation
==============

The official Oktoberfest documentation can be found at `https://oktoberfest.readthedocs.io <https://oktoberfest.readthedocs.io>`__.

How to cite
===========

.. include:: docs/reference.rst
   :start-line: 1
