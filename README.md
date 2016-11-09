# ScalarSolvers

[![Build Status](https://travis-ci.org/Goysa2/ScalarSolvers.jl.svg?branch=master)](https://travis-ci.org/Goysa2/ScalarSolvers.jl)

[![Coverage Status](https://coveralls.io/repos/Goysa2/ScalarSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Goysa2/ScalarSolvers.jl?branch=master)

[![codecov.io](http://codecov.io/github/Goysa2/ScalarSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/Goysa2/ScalarSolvers.jl?branch=master)

Collection of different algorithms to minimize problems like those in [ScalarOptimizationProblems](https://github.com/Goysa2/ScalarOptimizationProblems).

## How to use
With [NLPModels](https://github.com/JuliaSmoothOptimizers/NLPModels.jl), [ScalarOptimizationProblems](https://github.com/Goysa2/ScalarOptimizationProblems) and [Optimize](https://github.com/JuliaSmoothOptimizers/Optimize.jl) here is a simple example on how to use one of the algorithms:

`julia> model=MathProgNLPModel(AMPGO02())`

`julia> h=C1LineFunction(model,[0.0],[1.0])`

`julia> TR_Sec(h,2.7,7.5)`
