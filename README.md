# Fast Phase Factor Finding for QSP

This repository contains code for two efficient methods for finding phase factors in Quantum Signal Processing (QSP).

## Overview

The functions `HC` and `QSP_FFPI` implement two new methods for phase factor finding, as proposed in the paper [*Fast Phase Factor Finding for Quantum Signal Processing*](https://arxiv.org/abs/2410.06409). Additional auxiliary functions `half_chol`, `getUnitaryCoef`, and `getPQcoef` are included, along with selected functions from [QSPPACK](https://github.com/qsppack/QSPPACK) for preformance comparison purposes.

`HC` algorithm works for any *even* target polynomial $f$, even when $\max_{-1\le x \le 1}|f(x)|$ is close to 1. The code for odd polynomials is to be updated.

`QSP_FFPI` works for either even or odd target polynomial $f$, but may fail to converge when $\max_{-1\le x \le 1}|f(x)|$ is close to 1. In practice, $\max_{-1\le x \le 1}|f(x)|<\frac{1}{2}$ will make it work. When the degree `d` of the polynomial is greater than $10^5$, the default stop criteria $10^{-12}$ may be too strict as the machine error accumulates approximately linearly. A safe choice would be setting `opts.criteria = 1e-15*d`.

## Generating Figures

To replicate the figures presented in the paper, use the following scripts:

- `runtime_test_non_fully_coherent.m`: Generates Figure 1.
- `runtime_test_fully_coherent.m`: Generates Figure 2.
