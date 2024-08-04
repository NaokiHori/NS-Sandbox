# NS Sandbox

UNDER CONSTRUCTION

[![CI](https://github.com/NaokiHori/NS-Sandbox/actions/workflows/ci.yml/badge.svg)](https://github.com/NaokiHori/NS-Sandbox/actions/workflows/ci.yml)
[![License](https://img.shields.io/github/license/NaokiHori/NS-Sandbox)](https://opensource.org/license/MIT)
[![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/NS-Sandbox/main)](https://github.com/NaokiHori/NS-Sandbox/commits/main)

![thumbnail](https://github.com/NaokiHori/NS-Sandbox/blob/main/thumbnail.jpg)

An extremely simple two-dimensional incompressible Navier-Stokes solver designed to test various numerical algorithms.

## Concept

The primary goal of this library is simplicity, allowing numerical methods, especially for multi-phase flows, to be easily and quickly implemented and tested without being encumbered by sub-critical tasks. 
To achieve this, the following approaches are adopted:

- Euler-forward time-stepping
- Fully-explicit diffusive treatments
- No MPI parallelization

For the sake of transparency (and for fun), in-house Fourier transforms, linear matrix solvers, and matrix transpose routines are used, despite their sub-optimal performance.

## Dependencies

- [C Compiler](https://gcc.gnu.org)
- [GNU Make](https://www.gnu.org/software/make/)
- [OpenMP](https://www.openmp.org) (optional)

## Quick Start

1. Clone the repository:

    ```bash
    git clone https://github.com/NaokiHori/NS-Sandbox
    cd NS-Sandbox
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

Modify `include/domain.h` and re-build the source.

## Multi-thread Parallelization

Even in two-dimensional domains, some level of parallelization is necessary.
This project utilizes `OpenMP` for parallelization for convenience.

## Note

For simplicity, all flow fields have `NX + 2` by `NY + 2` elements, regardless of the type of arrays.
Here, `NX` and `NY` represent the degrees of freedom in the `x` and `y` directions, respectively.
This approach simplifies the incorporation of different boundary conditions (wall-bounded, periodic, inflow-outflow).

