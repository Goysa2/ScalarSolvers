export bissect_Cub
function bissect_Cub(h :: LineModel,
                    t1 :: Float64,
                    t2 :: Float64;
                    tol :: Float64 = 1e-7,
                    maxiter :: Int = 50,
                    verbose :: Bool = false)

  (tₐ, tᵦ) = trouve_intervalle(h, t1, 2.5)               

  γ = 0.8; t = tᵦ; tₚ = tₐ; tqnp = tₐ
  hₖ = 0; hkm1 = 0; gkm1 = 0; hₚ = 0
  iter = 0

  hₖ = obj(h, t); gₖ = grad(h, t)
  hkm1 = obj(h, tqnp); gkm1 = grad(h, tqnp)

  verbose &&
    @printf(" iter        tₚ        t         dN         ")
  verbose &&
    @printf("gₖ          gplus        \n")
  verbose &&
    @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
             iter, tₚ, t, 0.0, gₖ,  0.0)

  while ((abs(gₖ) > tol) && (iter < maxiter)) || (iter == 0)
    s = t - tqnp; y = gₖ - gkm1; α = -s
    z = gₖ + gkm1 + 3 * (hₖ - hkm1) / α
    discr = z^2 - gₖ * gkm1; denom = gₖ + gkm1 + 2 * z
    if (discr > 0.0) && (abs(denom) > eps(Float64))
      #si on peut on utilise l'interpolation cubique
      w = sqrt(discr)
      dN = -s * (gₖ + z + sign(α) * w) / (denom)
    else #on se rabat sur une étape de sécante
      dN = -gₖ * s / y
    end

    if ((tₚ - t) * dN > 0.0) & (dN / (tₚ - t) < γ)
      tplus = t + dN
      hplus = obj(h, tplus)
      gplus = grad(h, tplus)
      verbose && println("N")
    else
      tplus = (t + tₚ) / 2
      hplus = obj(h, tplus)
      gplus = grad(h, tplus)
      verbose && println("B")
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
  return (topt, iter)

end # function
