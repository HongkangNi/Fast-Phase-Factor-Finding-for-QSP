# Fast Phase Factor Finding for QSP

This repository contains code for two efficient methods for finding phase factors in Quantum Signal Processing (QSP).

## Overview

The functions `HC` and `QSP_FFPI` implement two new methods for phase factor finding, as proposed in the paper [*Fast Phase Factor Finding for Quantum Signal Processing*](https://arxiv.org/abs/2410.06409). Additional auxiliary functions `half_chol`, `getUnitaryCoef`, and `getPQcoef` are included, along with selected functions from [QSPPACK](https://github.com/qsppack/QSPPACK) for preformance comparison purposes.

## Generating Figures

To replicate the figures presented in the paper, use the following scripts:

- `runtime_test_non_fully_coherent.m`: Generates Figure 1.
- `runtime_test_fully_coherent.m`: Generates Figure 2.
