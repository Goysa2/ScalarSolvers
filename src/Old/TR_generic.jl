export TR_generic

# One dimensionnal trust region algorithm.
# We approximate our function h by a 2nd degree polynomial (Taylor) q.
# If q is a good approximation of h then we move towards the minimizer of q
# Otherwise we reduce our Trust region interval (Δ) until we have a good
# approximation.
# q(d) = h(t) + h'(t) * d + 0.5 * h''(t) * d²
#
# Specific key word arguments
# eps1 :: Float64. Below that treshold we don't move and reduce Δ.
#                  Treshold to determines if q is a bad approximation of h
# eps2 :: Float64. Over that treshold we augment the size of Δ
#                  Treshold to determines if q is a very good approximation of h
# red :: Float64. Rate at which we reduce Δ if we have a bad approximation
# aug :: Float64. Rate at which we augment Δ if we have a bad approximation
# Δ :: Float64. Size of the trust region at the begining of the algorithm.

function TR_generic(h :: AbstractNLPModel;
                    t₀ :: Float64 = h.meta.x0[1],
                    tol :: Float64 = 1e-7,
                    maxiter :: Int = 50,
                    verbose :: Bool = true,
                    eps1 :: Float64 = 0.2,
                    eps2 :: Float64 = 0.8,
                    red:: Float64 = 0.5,
                    aug :: Float64 = 2.0,
                    Δ :: Float64 = 1.0,
                    direction :: String = "Nwt")

    (length(h.meta.x0) > 1) && warn("Not a 1-D problem ")
    t = t₀; iter = 0; # Establish starting point and iteration counter.
    fₖ = obj(h, [t])[1]
    gₖ = grad(h, [t])[1]

    # H denotes the (approximation of the) second derivative
    if direction == "Sec" || direction == "SecA"
      H = 1.0
    else
      H = hess(h, [t])[1]
    end

    verbose &&
        @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter, t, gₖ[1], Δ)

    # We loop until we have a minimizer or we have reached the maximum number of
    # iterations.
    while ((abs(gₖ) > tol) && (iter < maxiter)) || (iter == 0)
        # @show gₖ
        # @show H
        dN = -gₖ / H;                 # (approximation of the) Newton direction
        # @show dN
        d = TR_step_computation(H, gₖ, dN, Δ) # determines the right direction

        # @show t                                      # depending on our trust region
        # @show d
        # @show t + d
        # @show [t + d]
        # Numerical reduction computation
        ftestTR = obj(h, [t + d])[1]        # Value of h and h' at t + d
        gtestTR = grad(h, [t + d])[1]
        # @show ftestTR
        # @show gtestTR

        # We check to see if we have a good approximation of h using the ratio
        # of the actual reduction and the predicted reduction.
        (pred, ared, ratio) =
            pred_ared_computation(gₖ, fₖ, H, d, ftestTR, gtestTR)
        # @show pred
        # @show ared
        # @show ratio
        # @show eps1
        # @show abs(d)
        # @show Δ
        if (ratio < eps1) && (abs(d) == Δ)
            # Bad approximation of h. Reduction of the size of the trust region
            Δ = red * Δ
            # printstyled("on a reduit delta \n", color = :red)
        else
            # Good approximation we move towards the minimizer of q
            (t, fₖ, gₖ, H) =
                step_computation(direction, h, t, d, gₖ, fₖ, ftestTR, gtestTR)
            # printstyled("après step computation \n", color = :red)
            # @show t
            # @show fₖ
            # @show gₖ
            # @show H
            if ratio > eps2
                # Very good approximation of h. We increase our trust region.
                Δ = aug * Δ
                # printstyled("on a augmente delta \n", color = :red)
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n",
                           iter, t, gₖ, Δ, pred, ared)
    end

    status = :NotSolved
    (abs(gₖ) < tol) && (status = :Optimal)
    (iter >= maxiter) && (status = :Tired)
    tired = iter > maxiter
    optimal = abs(gₖ) < tol
    return (t, fₖ, norm(gₖ, Inf), iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
