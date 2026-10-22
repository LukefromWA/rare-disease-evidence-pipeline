# Software environment

The recorded WSL library uses R 4.3.3. The [installed-package table](R-installed-packages.csv) preserves 677 package versions across computation, diagnostics, plotting and reporting. It is an installed-library inventory, not a minimal dependency list or a restoration-tested release lockfile.

The diagram generators use Python `diagrams` and Graphviz. Their pinned Python package versions and tool versions accompany the documentation sources in this directory.

The complete workflow environment also needs the selected Snakemake/Python environment, CmdStan and its compiler toolchain, and the reference/annotation tools invoked by the run. Do not substitute the package inventory for an explicit environment restoration specification.
