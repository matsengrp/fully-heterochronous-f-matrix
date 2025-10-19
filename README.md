# Fully Heterochronous F-Matrix

F-matrix representation of ranked tree shapes with Python and R implementations for sampling and enumeration.

## Overview

This repository contains code for working with F-matrix representations of ranked tree shapes, including both isochronous (all tips at the same time) and heterochronous (tips at different times) phylogenetic trees.

## Components

### Python (`python/`)
- **f_matrix.py**: Core functionality for F-matrix representation
  - Generate all F-matrices for given tree sizes
  - Convert between tree and matrix representations (F-matrix, D-matrix, E-matrix)
  - Sample from distributions on ranked tree shapes
  - Validate and visualize ranked tree shapes

### R (`Rscript/`)
- **enumerate.R / enumerate_clean.R**: Enumerate all isochronous trees
- **generate_fmatrices_optimized.R**: Optimized binary F-matrix generation
- **sampling.R**: Sampling utilities
- **simulation R code/**: Tree simulations
  - bernoulli_splitting.R
  - diag_top_down.R
  - sampling_coales.R

## Installation

```bash
# Create conda environment
mamba env create -f environment.yml
mamba activate farmrats

# Install package
pip install -e '.[dev]'
```

## Quick Start

```python
from python import f_matrix as f_mat

# Generate all F-matrices for 3-tip trees
matrices = f_mat.make_all_f_matrices(3)

# Convert tree to F-matrix
from ete3 import Tree
tree = Tree("((A:1,B:1):1,C:2);")
f_mat.tree_to_f_mat(tree)
```

## Testing

```bash
# Run tests
pytest tests/

# Run only F-matrix tests
pytest tests/test_f_matrix.py

# Run Python vs R comparison tests (requires R with ape package)
pytest tests/test_sampling.py -m r_comparison
```

## Repository

This code was extracted from the [farmrats](https://github.com/matsengrp/farmrats) repository, which contains deep learning models for ranked tree shapes.
