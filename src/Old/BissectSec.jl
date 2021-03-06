export bissect_sec
function bissect_sec(h :: AbstractNLPModel;
                     t1 :: Float64 = h.meta.x0[1],
                     t2 :: Float64 = 100.0,
                     tol :: Float64 = 1e-7,
                     maxiter :: Int = 50,
                     verbose :: Bool = false)

    (length(h.meta.x0) > 1) && warn("Not a 1-D problem ")
  (tₐ, tᵦ) = trouve_intervalle(h, t1, 2.5)
  γ = 0.8; t = tᵦ; tₚ = tₐ; tqnp = tₐ
  hₖ = 0; hkm1 = 0; gkm1 = 0; hplus = 0
  iter = 0

  gₖ = grad(h, [t])[1]
  gkm1 = grad(h, [tqnp])[1]

  verbose &&
    @printf(" iter        tₚ        t         dN         gₖ      ")
  verbose &&
    @printf("gplus        \n")
  verbose &&
    @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                iter, tₚ, t, 0.0, gₖ,  0.0)

  while ((abs(gₖ) > tol) && (iter < maxiter)) || (iter == 0)
    s = t - tqnp; y = gₖ - gkm1
    dN = -gₖ * s / y

    if ((tₚ - t) * dN > 0.0) && (dN / (tₚ - t) < γ)
      tplus = t + dN
      #hplus = obj(h, tplus)
      gplus = grad(h, [tplus])[1]
      verbose && print_with_color(:green, "N")
    else
      tplus = (t + tₚ) / 2
      #hplus = obj(h,  tplus)
      gplus = grad(h, [tplus])[1]
      verbose && print_with_color(:green, "B")
    end

    if t > tₚ
      if gplus < 0.0
        tₚ = t; tqnp = t; t = tplus
      else
        tqnp = t; t = tplus
      end
    else
      if gplus > 0.0
        tₚ = t; tqnp = t; t = tplus
      else
        tqnp = t; t = tplus
      end
    end

    #mise à jour des valeurs
    hkm1 = hₖ; gkm1 = gₖ; hₖ = hplus; gₖ = gplus
    iter = iter + 1
    verbose &&
        @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                iter, tₚ, t, dN, gₖ, gplus)
  end
    topt = t
    status = :NotSolved
    (abs(gₖ) < tol) && (status = :Optimal)
    (iter >= maxiter) && (status = :Tired)
    tired = iter > maxiter
    optimal = abs(gₖ) < tol
    return (topt, obj(h, [topt])[1], norm(gₖ, Inf), iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
