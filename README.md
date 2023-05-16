# FewChallenger

How blazingly fast is FEW?

In this project, we try to find speed-up opportunities in the low-level C layer of [FEW](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms).
We focus on the amplitude computation for a first hands on session.
The overall memory layout is challenged in [MemoryChallenger](MemoryChallenger/MemoryChallenger)
and the spline interpolation implementation in [SplineChallenger](SplineChallenger/SplineChallenger).

A benchmark is implemented as executable [Challenge](Challenger/src/program/Challenge.cpp).
Preliminary results show a ~20% speed-up with memory optimization only (~5s vs. ~4s for 1000 trajectory points),
and a factor ~6 improvement with a complete reorganization of the computation (~0.8s).

# FEW

The amplitude `A_lmn(p, e)` is a complex-valued function.
It is known for some values (the knots) of its parameters and has to be interpolated over a vector of `(p, e)`'s.
The so-called modes are such that:

* `l` ranges from 2 to `lmax`;
* `m` ranges from 0 to `l`;
* `n` ranges from `-nmax` to `nmax`.

In the original code, one pair of real interpolants (of type `gsl_spline`) is instantiated for each mode, to account for the real and imaginary parts.
They are organized in some chain of pointers `gsl_spline***`.

Each pair of interpolants is then applied to each `(p, e)` to get the amplitude real and imaginary parts.
For each mode, the interpolant is called on the same positions, i.e. there is a single vector of `(p, e)`'s for all the `(l, m, n)`'s.

The default benchmark uses the following sizing parameters:

* `lmax` = 10;
* `nmax` = 30;
* There are 33 x 50 knots;
* There are 1000 interpolation coordinates.

# MemoryChallenger

The memory challenge consists in rearranging the interpolants in contiguous memory similar to `Interpolant[]`,
with some indexing scheme which does the mapping between 1D index in the array and `(l, m, n)`.

For simplicity, this is done for now in some home-made ndarray ([`Linx::Raster`](https://github.com/kabasset/Linx)).
It has a box domain, such that half the storage is wasted (where `m > l`).
Still, data is better packed and cache friendlier.

Without changing the actual computation, this refactoring of the memory layout saves around 20% of computation time.

# SplineChallenger

Instead of using a pair of real interpolants for the real and imaginary parts of the amplitude, one can directly create complex interpolants.
Given that most of the spline coefficients are real anyway, this should almost divide the computation time by 2.

Additionally, given that the grid of positions to interpolate on does *not* depend on the mode, many spline coefficients can be computed intependently of `(l, m, n)`.
Let us name `U, V_lmn` the knot positions and values, and `X, Y_lmn` the interpolated positions and values.
Instead of computing one interpolant for each `(l, m, n)` and applying it to the `X`'s, we propose to compute a single interpolant -- or, rather, *resampler* -- for the `U`'s and `X`'s.
A lot of spline coefficients can be computed already.
Moreover, given that `X` is known early, only the knots which neighbor `X` values have to be computed, which is not possible with classical spline interpolation.
Finally, only the coefficients which depend on `V_lmn` around `X` must be computed to get `Y_lmn` for each mode.
Iteratively loading `V_lmn` and computing `Y_lmn` for each mode allows instanciating a single interpolant.

Unfortunately, we did not find a suitable library and therefore had to implement one: [Splider](https://github.com/kabasset/Splider).

With respect to the original code, the use of Splider instead of GSL divides computation time by around 6.
