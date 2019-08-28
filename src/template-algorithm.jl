# Generic template for 1-D optimization.
# Used only to explain the parameters. This function does absolutely nothing and
# shouldn't be used. Ever.
#
# We are minimizing h(t): R → R.
#
# Inputs (common to all algorithm, specifics inputs/hyperparameters are
# explained in each function):
# h       :: Line Model.  Line model is a tool used for Line Search. We can adapt it
# 	   		 to one dimensionnal optimization. See the README to know how to create
#      		 a Line model
# nlpstop :: A Stopping object dealing with the stopping criterion of the
#    		 algorithm. 
# verbose :: Bool. If true, print iteration of linesearch

function one_dimensional_template(h       :: AbstractNLPModel,
                         		  nlpstop :: AbstractStopping;
			                      verbose :: Bool = true,
								  kwargs...)
	# Core of the function

	# Outputs:
	# optimal :: Status if the problem has been solved
	# nlpstop :: Information about the output

	optimal = false;
	return optimal, nlpstop
end
