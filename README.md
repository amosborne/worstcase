# worstcase

`pip install worstcase`

Worst case analysis and sensitivity studies using Extreme Value and/or Monte Carlo methods.

EDN blogger Charles Hymowitz wrote an insightful series on worst case circuit analysis, read it [here](https://www.edn.com/the-worst-case/).

This lightweight package allows the user to specify varying parameters and execute worst case analysis or sensitivity studies by Extreme Value and/or Monte Carlo methods over arbitrary single valued functions. Parameters are assumed to be uniform by default; the user may generate their own custom distributions if desired.

Please see the [example](./examples/voltage_divider.ipynb) for usage.