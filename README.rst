Oktoberfest: Rescoring, Collision Energy Calibration and Spectral Library Generation for Proteomics
===================================================================================================

|PyPI| |Python Version| |License| |Read the Docs| |Build| |Tests| |Codecov| |pre-commit| |Black|

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
.. |Build| image:: https://github.com/wilhelm-lab/oktoberfest/workflows/Build%20oktoberfest%20Package/badge.svg
   :target: https://github.com/wilhelm-lab/oktoberfest/actions?workflow=Package
   :alt: Build Package Status
.. |Tests| image:: https://github.com/wilhelm-lab/oktoberfest/workflows/Run%20oktoberfest%20Tests/badge.svg
   :target: https://github.com/wilhelm-lab/oktoberfest/actions?workflow=Tests
   :alt: Run Tests Status
.. |Codecov| image:: https://codecov.io/gh/wilhelm-lab/oktoberfest/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/wilhelm-lab/oktoberfest
   :alt: Codecov
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
   :alt: pre-commit
.. |Black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
   :alt: Black

Oktoberfest is a python tool for collision energy calibration, rescoring search results and generating spectral libraries for proteomics research within the Prosit ecosystem. It offers an end to end pipeline that takes search results, predicts peptide properties using koina, plots summaries and quality control figures and performs FDR estimation with either mokapot or percolator.
