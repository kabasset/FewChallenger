# FewChallenger

How blazingly fast is FEW?

# Purpose

In this project, we try to find speed-up opportunities in the low-level C layer of [FEW](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms).
We focus on the amplitude computation for a first hands on session.
It relies on some spline interpolation to reduce complexity drastically.

The amplitude `A_lmn(p, e)` is a complex-valued function parametrized by some mode indices `l`, `m` and `n`.
It is known for some values (the knots) of its parameters and has to be interpolated over a vector of `(p, e)`'s named trajectory.
The modes are such that:

* `l` ranges from 2 to `lmax`;
* `m` ranges from 0 to `l`;
* `n` ranges from `-nmax` to `nmax`.

The overall memory layout is challenged in [MemoryChallenger](MemoryChallenger/MemoryChallenger),
where we target better data locallity.
The spline interpolation implementation itself is challenged in [SplineChallenger](SplineChallenger/SplineChallenger),
where we cache many interediate results which can be reused for every mode.
More details are given in the next sections.

# Benchmark

A benchmark is implemented as executable [Challenge](Challenger/src/program/Challenge.cpp).
By default, the benchmark uses the following sizing parameters:

* `lmax` = 10, `nmax` = 30, which yields 3843 modes;
* There are 33 x 50 = 1650 knots;
* There are 1000 trajectory points.

Here are preliminary results:

| Implementation      | Runtime | Speedup |
| ------------------- | ------- | ------- |
| Original FEW        |    4.4s |     1.0 |
| Memory optimization |    3.9s |     1.1 |
| Spline caching      |   0.09s |  **50** |

# FEW

In the original code, the interpolation relies on the GSL.
One pair of real interpolants (of type `gsl_spline`) is instantiated for each mode, to account for the real and imaginary parts.
They are organized in some chain of pointers `gsl_spline***`.

Each pair of interpolants is then applied to each `(p, e)` to get the amplitude real and imaginary parts.
For each mode, the interpolant is called on the same positions, i.e. there is a single vector of `(p, e)`'s for all the `(l, m, n)`'s.

# MemoryChallenger

The memory challenge consists in rearranging the interpolants in contiguous memory similar to `Interpolant[]` instead of FEW's `Interpolant***`,
with some indexing scheme which does the mapping between 1D index in the array and `(l, m, n)`.

For simplicity, this is done for now in some home-made ndarray ([`Linx::Raster`](https://github.com/kabasset/Linx)).
It has a box domain, such that half the storage is wasted (where `m > l`).
Still, data is better packed and CPU cache-friendlier.

# SplineChallenger

Instead of using a pair of real interpolants for the real and imaginary parts of the amplitude, one can directly create complex interpolants
(which the GSL does not support).
Given that most of the spline coefficients are real anyway, this should almost divide the computation time by 2.

Additionally, given that the grid of positions to interpolate on does *not* depend on the mode, many spline coefficients can be computed intependently of `(l, m, n)`.
Let us name `U, V_lmn` the knot positions and values, and `X, Y_lmn` the interpolated positions and values.
Instead of computing one interpolant for each `(l, m, n)` and applying it to `X`, we propose to compute a single parametric interpolant from `U` to `X`, and to apply it to the `V_lmn`'s.
A lot of spline coefficients can be computed already.
Moreover, given that `X` is known early, only the knots which neighbor `X` values have to be computed, which is not possible with classical spline interpolation.
Finally, only the coefficients which depend on `V_lmn` around `X` must be computed to get `Y_lmn` for each mode.
Iteratively loading `V_lmn` and computing `Y_lmn` for each mode allows instanciating a single interpolant.

Unfortunately, we did not find a suitable library and therefore had to implement one: [Splider](https://github.com/kabasset/Splider).
Many optimization are still to be implemented in Splider, but the speedup is already spectacular.
