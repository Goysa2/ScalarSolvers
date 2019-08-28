# ScalarSolvers

[![Build Status](https://travis-ci.org/Goysa2/ScalarSolvers.jl.svg?branch=master)](https://travis-ci.org/Goysa2/ScalarSolvers.jl)

[![Coverage Status](https://coveralls.io/repos/Goysa2/ScalarSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Goysa2/ScalarSolvers.jl?branch=master)

[![codecov.io](http://codecov.io/github/Goysa2/ScalarSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/Goysa2/ScalarSolvers.jl?branch=master)

Collection of different algorithms to minimize problems like the one dimensionnal problems in [OptimizationProblems](https://github.com/JuliaSmoothOptimizers/OptimizationProblems.jl).

## How to use
With [NLPModels](https://github.com/JuliaSmoothOptimizers/NLPModels.jl), [OptimizationProblems](https://github.com/JuliaSmoothOptimizers/OptimizationProblems.jl), [State](https://github.com/Goysa2/State.jl) and [Stopping](https://github.com/Goysa2/Stopping.jl). Here is a simple example:
```
nlp = MathProgNLPModel(AMPGO02());
nlp_at_x = NLPAtX([3.7])
nlp_stop = NLPStopping(nlp, Stopping.unconstrained, nlp_at_x)
optimal, stop = TR_Nwt(nlp, nlp_stop, verbose = true)
```
The output should be:

```
 iter  t         gₖ          Δ        pred         ared
    0 2.70e+00  -3.94e+00  1.00e+00
    1 3.70e+00  2.40e+00  1.00e+00 -6.44e+00 -1.60e+00
    2 2.93e+00  -4.13e+00  1.00e+00 -9.27e-01 6.51e-01
    3 2.93e+00  -4.13e+00  5.00e-01 -2.44e+00 -9.85e-02
    4 3.43e+00  4.16e-01  5.00e-01 -1.64e+00 -1.08e+00
    5 3.39e+00  -1.02e-02  1.00e+00 -8.33e-03 -8.18e-03
    6 3.39e+00  -4.32e-06  2.00e+00 -4.83e-06 -4.83e-06
    7 3.39e+00  -7.75e-13  4.00e+00 -8.56e-13 -8.56e-13

```

## Installing
`julia> add https://github.com/Goysa2/ScalarSolvers`

## Line Search applications
All the algorithms presented in this packages are adapted for Line Search in
the package [LineSearch](https://github.com/Goysa2/LineSearch).
It is possible to directly use these algorithm for line search using an appropriate stopping criterion with a Stopping object.

## References
All the algorithms can be found here:

J. Nocedal and S.Wright [Numerical Optimization](http://www.bioinfo.org.cn/~wangchao/maa/Numerical_Optimization.pdf)

J.P Dusseault Univariate diffentiable optiization algorithms and LineSearch
computation
