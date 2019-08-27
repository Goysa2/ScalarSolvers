export ARC_generic_stop

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

function ARC_generic_stop(h         :: AbstractNLPModel,
                          nlpstop   :: AbstractStopping;
                          t₀        :: Float64 = h.meta.x0[1],
                          tol       :: Float64 = 1e-7,
                          maxiter   :: Int = 50,
                          verbose   :: Bool = true,
                          eps1      :: Float64 = 0.25,
                          eps2      :: Float64 = 0.75,
                          red       :: Float64 = 0.2,
                          aug       :: Float64 = 5.0,
                          Δ         :: Float64 = 0.5,
                          direction :: Symbol = :Nwt)

    Δ = [Δ]
    t = h.meta.x0; iter = 0;        # We establish our starting point t
    fₖ = obj(h, t);
    f₀ = copy(fₖ)
    gₖ = grad(h, t);
    g₀ = copy(gₖ)


    # H denotes the (approximation of the) second derivative
    if direction == :Nwt
      H = hess(h, t)
  elseif direction == :Sec || direction == :SecA
      H = [1.0]
    end

    OK = update_and_start!(nlpstop, x = t, fx = fₖ, gx = gₖ, g0 = g₀, Hx = H)

    verbose &&
        @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose &&
        @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter, t[1], gₖ[1], Δ[1])
    # We loop until we have a minimizer or we have reached the maximum number of
    # iterations.
    while !OK

        d = ARC_step_computation_stop(H, gₖ, Δ) # We find the direction in which we
                                                # move
        # @show t
        # @show d
        # Numerical reduction computation
        ftestTR = obj(h, t + d)  # Value of h and h' at t + d
        gtestTR = grad(h, t + d)

        # We check to see if we have a good approximation of h using the ratio
        # of the actual reduction and the predicted reduction.
        (pred, ared, ratio) =
            pred_ared_computation_stop(gₖ, fₖ, H, d, ftestTR, gtestTR)

        if (ratio .< eps1)
            # Bad approximation of h. We make the cubic term more prevalant
            Δ = red * Δ
        else
            # Good approximation we move towards the minimizer of c
            (t, fₖ, gₖ, H) =
                step_computation_stop(direction, h, t, d, gₖ, fₖ, ftestTR, gtestTR)

            if ratio > eps2
                # Very good approximation of h. We reduce the importance of the
                # cubic term.
                Δ = aug * Δ
            end
            OK = update_and_stop!(nlpstop, x = t, fx = fₖ, gx = gₖ, Hx = H)
        end

        iter += 1
        # @show iter
        # @show t[1]
        # @show gₖ[1]
        # @show Δ[1]
        # @show pred[1]
        # @show ared[1]
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n",
                            iter, t[1], gₖ[1], Δ[1], pred[1], ared[1])
    end

    optimal = OK
    return optimal, nlpstop
end
