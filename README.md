# Very Simple NS Solver

[![CI](https://github.com/NaokiHori/VerySimpleNSSolver/actions/workflows/ci.yml/badge.svg)](https://github.com/NaokiHori/VerySimpleNSSolver/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/NaokiHori/VerySimpleNSSolver)](https://opensource.org/license/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/VerySimpleNSSolver/main)](https://github.com/NaokiHori/VerySimpleNSSolver/commits/main)

[![cover](https://github.com/NaokiHori/VerySimpleNSSolver/blob/main/cover.jpg)](https://youtu.be/MdViM5HzCDs)

An extremely simple two-dimensional incompressible Navier-Stokes solver designed to test various numerical algorithms.

## Concept

The primary goal of this library is simplicity, allowing numerical methods, especially for multi-phase flows, to be easily and quickly implemented and tested without being encumbered by sub-critical tasks. 
To achieve this, the following approaches are adopted:

- Euler-forward time-stepping
- Fully-explicit diffusive treatments
- No multi-process (e.g., `MPI`) parallelization

For the sake of transparency (and for fun), in-house Fourier transforms, linear matrix solvers, and matrix transpose routines are used, despite their sub-optimal performance.

## Dependencies

- [C Compiler](https://gcc.gnu.org)
- [GNU Make](https://www.gnu.org/software/make/)
- [OpenMP](https://www.openmp.org) (optional)

## Quick Start

1. Clone the repository:

    ```bash
    git clone https://github.com/NaokiHori/VerySimpleNSSolver
    cd VerySimpleNSSolver
    ```

2. Build the project:

    ```bash
    make all
    ```

    If your platform supports `OpenMP`:

    ```bash
    make ARG_CFLAG="-fopenmp" all
    ```

    Alternatively, you can directly modify the `Makefile`.

3. Run the solver:

    ```bash
    make output
    ./a.out
    ```

    If your platform supports `OpenMP`:

    ```bash
    make output
    OMP_NUM_THREADS=4 ./a.out
    ```

    (This example uses four threads.)

## Domain Size

The domain sizes (spatial resolutions and lengths) are defined in `src/domain.c`.
Modify the corresponding parameter and re-build the source.

## Multi-thread Parallelization

Even in two-dimensional domains, some level of parallelization is necessary.
This project utilizes `OpenMP` for parallelization for convenience.

## Note

For simplicity, all flow fields have `domain->nx + 2` by `domain->ny + 2` elements, regardless of the type of arrays.
Here, `domain->nx` and `domain->ny` represent the degrees of freedom in the `x` and `y` directions, respectively.
This approach simplifies the incorporation of different boundary conditions (wall-bounded, periodic, inflow-outflow).

