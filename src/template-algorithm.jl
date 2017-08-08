# Generic template for 1-D optimization.
# Used only to explain the parameters. This function does absolutely nothing and
# shouldn't be used. Ever.
#
# We are minimizing h(t): R → R.
#
# Inputs (common to all algorithm, specifics inputs/hyperparameters are
# explained in each function):
# h :: Line Model.  Line model is a tool used for Line Search. We can adapt it
# 	   to one dimensionnal optimization. See the README to know how to create
#      a Line model
# t₀ :: Float64. Starting point of the algorithm.
# tₘ :: Float64. Maximum bound. If a local minimizer exists it will be found
#       in (t₀, tₘ)
# verbose :: Bool. If true, print iteration of linesearch
# maxiter :: Int64. Maximum number of iteration of the algorithm.
# tol :: Float64. Precision we set to determine a local minimizer

function one_dimensional_template(h :: LineModel,
                     			  t₀ :: Float64,
 	   		                      tₘ :: Float64;
	   		                      tol :: Float64 = 1e-7,
			                      maxiter :: Int = 50,
			                      verbose :: Bool = true,
								  kwargs...)
	# Core of the function

	# Outputs:
	# topt :: Float64. Our t* such as h(t*) minimize h.
	# iter :: Number of iteration of the algorithm

	topt = NaN; iter = NaN;
	return (topt, iter)
end
