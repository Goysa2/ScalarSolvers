export ARC_generic

# One dimensionnal ARC algorithm.
# We approximate our function h by a 2nd degree polynomial (Taylor) q enriched
# with a cubic term.
# If c is a good approximation of h then we move towards the minimizer of c
# Otherwise we reduce Δ in order to increase the weight of the cubic term.
# c(d) = h(t) + h'(t) * d + 0.5 * h''(t) * d² + (1/(3*Δ))*|d|³
#
# Specific key word arguments
# eps1 :: Float64. Below that treshold we don't move and reduce Δ.
#                  Treshold to determines if c is a bad approximation of h
# eps2 :: Float64. Over that treshold we augment the size of Δ
#                  Treshold to determines if c is a very good approximation of h
# red :: Float64. Rate at which we reduce Δ if we have a bad approximation
# aug :: Float64. Rate at which we augment Δ if we have a bad approximation
# Δ :: Float64. Size of the trust region at the begining of the algorithm.

function ARC_generic(h :: LineModel,
                     t₀ :: Float64,
                     tₘ :: Float64;
                     tol :: Float64 = 1e-7,
                     maxiter :: Int = 50,
                     verbose :: Bool = true,
                     eps1 :: Float64 = 0.25,
                     eps2 :: Float64 = 0.75,
                     red :: Float64 = 0.2,
                     aug :: Float64 = 5.0,
                     Δ :: Float64 = 0.5,
                     direction :: String = "Nwt")

    t = t₀; iter = 0;                # We establish our starting point t
    fₖ = obj(h, t); gₖ = grad(h, t); # And h(t) and h''(t)

    # H denotes the (approximation of the) second derivative
    if direction == "Nwt"
      H = hess(h, t)
    elseif direction == "Sec" || direction == "SecA"
      H = 1.0
    end

    #q(d) = fₖ + gₖ*d + 0.5*secₖ*d^2 + (1/3*(Δ))*abs(d)^3

    verbose &&
        @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose &&
        @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter, t, gₖ, Δ)
    # We loop until we have a minimizer or we have reached the maximum number of
    # iterations.
    while ((abs(gₖ) > tol) & (iter < maxiter)) | (iter == 0)

        d = ARC_step_computation(H, gₖ, Δ) # We find the direction in which we
                                           # move
        # Numerical reduction computation
        ftestTR = obj(h, t + d)  # Value of h and h' at t + d
        gtestTR = grad(h, t + d)

        # We check to see if we have a good approximation of h using the ratio
        # of the actual reduction and the predicted reduction.
        (pred, ared, ratio) =
            pred_ared_computation(gₖ, fₖ, H, d, ftestTR, gtestTR)

        if (ratio < eps1)
            # Bad approximation of h. We make the cubic term more prevalant
            Δ = red * Δ
        else
            # Good approximation we move towards the minimizer of c
            (t, fₖ, gₖ, H) =
                step_computation(direction, h, t, d, gₖ, fₖ, ftestTR, gtestTR)

            if ratio > eps2
                # Very good approximation of h. We reduce the importance of the
                # cubic term. 
                Δ = aug * Δ
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n",
                            iter, t, gₖ, Δ, pred, ared)
    end

    return (t, iter)
end
