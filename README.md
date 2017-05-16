# ScalarSolvers

[![Build Status](https://travis-ci.org/Goysa2/ScalarSolvers.jl.svg?branch=master)](https://travis-ci.org/Goysa2/ScalarSolvers.jl)

[![Coverage Status](https://coveralls.io/repos/Goysa2/ScalarSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Goysa2/ScalarSolvers.jl?branch=master)

[![codecov.io](http://codecov.io/github/Goysa2/ScalarSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/Goysa2/ScalarSolvers.jl?branch=master)

Collection of different algorithms to minimize problems like those in [ScalarOptimizationProblems](https://github.com/Goysa2/ScalarOptimizationProblems).

## How to use
With [NLPModels](https://github.com/JuliaSmoothOptimizers/NLPModels.jl), [ScalarOptimizationProblems](https://github.com/Goysa2/ScalarOptimizationProblems) and [Optimize](https://github.com/JuliaSmoothOptimizers/Optimize.jl) here is a simple example on how to use one of the algorithms:

`julia> model=MathProgNLPModel(AMPGO02());`

`julia> h=C2LineFunction(model,[0.0],[1.0]);`

`julia> new_TR_generic(h,2.7,7.5)`

By default the Newton direction will be used.

The output should be:


`iter  t         gₖ          Δ        pred         ared`

`   0 2.70e+00  -3.94e+00  1.00e+00 `

`   1 3.70e+00  2.40e+00  1.00e+00 -6.44e+00 -1.60e+00`

`   2 3.70e+00  2.40e+00  5.00e-01 -9.27e-01 6.51e-01`

`   3 3.20e+00  -2.07e+00  5.00e-01 -8.11e-01 -2.44e-01`

`   4 3.40e+00  9.73e-02  1.00e+00 -2.04e-01 -1.95e-01`

`   5 3.39e+00  -4.28e-04  2.00e+00 -4.39e-04 -4.38e-04`

`   6 3.39e+00  -7.60e-09  4.00e+00 -8.43e-09 -8.43e-09`

`(3.3872517177456873,6)`

## Installing
`julia> Pkg.clone("https://github.com/Goysa2/ScalarSolvers.git")`

## Sources
All the algorithms can be found here:

J. Nocedal and S.Wright [Numerical Optimization](http://www.bioinfo.org.cn/~wangchao/maa/Numerical_Optimization.pdf)

J.P Dusseault Univariate diffentiable optiization algorithms and LineSearch computation
