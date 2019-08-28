# This is a template for the bisection algorithm since they all follow the same
# pattern. This function doesn't do anything.
#
# Same inputs as other algorithms

function template_bisection(h :: LineModel,
                    		t1 :: Float64,
		                    t2 :: Float64;
		                    tol :: Float64 = 1e-7,
		                    maxiter :: Int = 50,
		                    verbose :: Bool = false,
							kwargs...)
	# Find an interval containing a local minimizer of h
	(tₐ, tᵦ) = trouve_intervalle(h, t1, 2.5)

	# Until we reach a minimizer or cut the interval in half or move towards
	# the Newton direction.

	return (topt, iter)
end
