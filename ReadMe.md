[![PyPI](https://img.shields.io/pypi/v/oktoberfest.svg)](https://pypi.org/project/oktoberfest/)
[![Python Version](https://img.shields.io/pypi/pyversions/oktoberfest)](https://pypi.org/project/oktoberfest)
[![License](https://img.shields.io/github/license/wilhelm-lab/oktoberfest)](https://opensource.org/licenses/MIT)
[![Read the Docs](https://img.shields.io/readthedocs/oktoberfest/latest.svg?label=Read%20the%20Docs)](https://oktoberfest.readthedocs.io/)
[![Build](https://github.com/wilhelm-lab/oktoberfest/workflows/Build%20oktoberfest%20Package/badge.svg)](https://github.com/wilhelm-lab/oktoberfest/actions?workflow=Package)
[![Tests](https://github.com/wilhelm-lab/oktoberfest/workflows/Run%20oktoberfest%20Tests/badge.svg)](https://github.com/wilhelm-lab/oktoberfest/actions?workflow=Tests)
[![Codecov](https://codecov.io/gh/wilhelm-lab/oktoberfest/branch/main/graph/badge.svg)](https://codecov.io/gh/wilhelm-lab/oktoberfest)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# Oktoberfest: Rescoring, Collision Energy Calibration and Spectral Library Generation for Proteomics

Oktoberfest is a python tool for collision energy calibration, rescoring search results and generating spectral libraries for proteomics research within the Prosit ecosystem. It offers an end to end pipeline that takes search results, predicts peptide properties using koina, plots summaries and quality control figures and performs FDR estimation with either mokapot or percolator.

## Documentation

The official Oktoberfest documentation can be found at https://oktoberfest.readthedocs.io.

## How to cite

### Please always cite the main publication

[Oktoberfest] Picciani M, Gabriel W, Giurcoiu VG et al. (2023), _Oktoberfest: Open-source spectral library generation and rescoring pipeline based on Prosit_, [Proteomics](https://doi.org/10.1002/pmic.202300112)

### Please cite the Prosit model(s) you are using

When using Prosit

[Prosit] Gessulat S, Schmidt T, Zolg DP et al. (2019), _PROSIT: Proteome-wide prediction of peptide tandem mass spectra by deep learning_, [Nature Methods](https://doi.org/10.1038/s41592-019-0426-7)

[Prosit-HLA] Wilhelm M, Zolg DP, Graber M et al. (2021), _Deep learning boosts sensitivity of mass spectrometry-based immunopeptidomics_, [Nature Communications](https://doi.org/10.1038/s41467-021-23713-9)

When using Prosit-TMT

[Prosit-TMT] Gabriel W, The M, Zolg D et al. (2022), _Prosit-TMT: Deep Learning Boosts Identification of TMT-Labeled Peptides_, [Analytical Chemistry](https://doi.org/10.1021/acs.analchem.1c05435)

When using Prosit-timsTOF

[Prosit-timsTOF] Adams C, Gabriel W, Laukens K et al. (2024), _Fragment ion intensity prediction improves the identification rate of non-tryptic peptides in timsTOF_, [Nature Communications](https://doi.org/10.1038/s41467-024-48322-0)

When using Prosit-XL

[Prosit-XL] Kalhor M, Saylan C, Picciani M et al. (2024), _Prosit-XL: enhanced cross-linked peptide identification by accurate fragment intensity prediction to study protein-protein interactions and protein structures_, [BioRxiv](https://www.biorxiv.org/content/10.1101/2024.12.15.627797v1)

### Please cite when using protein grouping and quantification

[PickedGroupFDR] The M, Samaras P, Kuster B, Wilhelm, M. (2022), _Reanalysis of ProteomicsDB using an accurate, sensitive, and scalable false discovery rate estimation approach for protein groups_, [Molecular & Cellular Proteomics](https://doi-org.org/10.1016/j.mcpro.2022.100437)
