# FewChallenger

How blazingly fast is FEW?

In this project, we try to find speed-up opportunities in the low-level C layer of [FEW](https://github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms).
We focus on the amplitude computation for a first hands on session.
The overall memory layout is challenged in [MemoryChallenger](MemoryChallenger/MemoryChallenger)
and the spline interpolation implementation in [SplineChallenger](SplineChallenger/SplineChallenger).

A benchmark is implemented as executable [Challenge](Challenger/src/program/Challenge.cpp).
Preliminary results show a ~20% speed-up with memory optimization (~5s vs. ~4s for 1000 trajectory points).
