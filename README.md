# worstcase

## Overview

`pip install worstcase`

Worst case analysis and sensitivity studies using Extreme Value and/or Monte Carlo methods.

This package coexists alongside far more capable uncertainty analysis and error propagation packages such as [uncertainties](https://pypi.org/project/uncertainties/) (first-order propagation), [soerp](https://pypi.org/project/soerp/) (second-order propagation), and [mcerp](https://pypi.org/project/mcerp/) (Monte Carlo propagation).

This package is designed for engineering applications where worst case analysis computations are often done using the Extreme Value method over single-valued functions while falling back to the Monte Carlo method when the worst case is known to not exist at the extremes. The Extreme Value method is implemented as a brute-force search; the Monte Carlo method is implemented with Latin Hypercube Sampling over a uniform distribution.

## Usage

See the example usage (here)[https://github.com/amosborne/worstcase/blob/master/examples/readme.ipynb].
