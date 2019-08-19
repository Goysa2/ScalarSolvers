export TR_generic_Stop

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

function TR_generic_Stop(h :: AbstractNLPModel,
                         nlpstop :: AbstractStopping;
                         verbose :: Bool = true,
                         eps1 :: Float64 = 0.2,
                         eps2 :: Float64 = 0.8,
                         red:: Float64 = 0.5,
                         aug :: Float64 = 2.0,
                         Δ :: Float64 = 1.0,
                         direction :: Symbol = :Nwt)

    (length(h.meta.x0) > 1) && warn("Not a 1-D problem ")
    t = h.meta.x0; iter = 0; # Establish starting point and iteration counter.
    fₖ = obj(h, t)
    f₀ = copy(fₖ)
    gₖ = grad(h, t)
    g₀ = copy(gₖ)
    printstyled("on a f et g \n", color = :green)

    # H denotes the (approximation of the) second derivative
    if direction == "Sec" || direction == "SecA"
      H = [1.0]
    else
      H = hess(h, t)
    end

    OK = update_and_start!(nlpstop, x = t, fx = fₖ, gx = gₖ, g0 = g₀, Hx = H)
    @show OK


    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter, t[1], gₖ[1], Δ)

    # We loop until we have a minimizer or we have reached the maximum number of
    # iterations.
    while !OK

        dN = -gₖ ./ H;                 # (approximation of the) Newton direction
        # @show dN
        d = TR_step_computation_stop(H, gₖ, dN, Δ) # determines the right direction
                                              # depending on our trust region
        # printstyled("on a d = $d \n", color = :red)
        # @show t
        # @show d
        # Numerical reduction computation
        ftestTR = obj(h, t + d)        # Value of h and h' at t + d
        gtestTR = grad(h, t + d)
        # @show ftestTR
        # @show gtestTR

        # We check to see if we have a good approximation of h using the ratio
        # of the actual reduction and the predicted reduction.
        (pred, ared, ratio) =
            pred_ared_computation_stop(gₖ, fₖ, H, d, ftestTR, gtestTR)
        # printstyled("une fois calculé \n", color = :bold)
        # @show pred
        # @show ared
        # @show ratio

        if (ratio .< eps1) && (abs.(d) == Δ)
            # Bad approximation of h. Reduction of the size of the trust region
            Δ = red * Δ
        else
            # Good approximation we move towards the minimizer of q
            (t, fₖ, gₖ, H) =
                step_computation_stop(direction, h, t, d, gₖ, fₖ, ftestTR, gtestTR)

            if ratio > eps2
                # Very good approximation of h. We increase our trust region.
                Δ = aug * Δ
            end
        end

        OK = update_and_stop!(nlpstop, x = t, fx = fₖ, gx = gₖ, Hx = H)

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n",
                           iter, t[1], gₖ[1], Δ, pred, ared)
    end

    status = :NotSolved
    OK && (status = :Optimal)
    # (iter >= maxiter) && (status = :Tired)
    # tired = iter > maxiter
    optimal = OK #abs(gₖ) < tol
    return (t, fₖ, norm(gₖ, Inf), iter, optimal, false, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
