export bissect
function bissect_Cub(h :: AbstractNLPModel,
                     nlpstop :: AbstractStopping;
                     t1 = h.meta.x0, t2 = 100.0,
                     verbose = false)


  (tₐ, tᵦ) = trouve_intervalle(h, t1, 2.5)

  γ = 0.8; t = tᵦ; tₚ = tₐ; tqnp = tₐ
  hₖ = [0.0]; hkm1 = [0.0]; gkm1 = [0.0]; hₚ = [0.0]
  iter = 0

  hₖ   = obj(h, t);    gₖ   = grad(h, t)
  hkm1 = obj(h, tqnp); gkm1 = grad(h, tqnp)

  OK = update_and_start!(nlpstop, x = t, gx = gₖ, g0 = copy(gₖ))

  verbose &&
    @printf(" iter        tₚ        t         dN         ")
  verbose &&
    @printf("gₖ          gplus        \n")
  verbose &&
    @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
             iter, tₚ[1], t[1], 0.0, gₖ[1],  0.0)

  while !OK
    s = t .- tqnp; y = gₖ .- gkm1; α = -s
    z = gₖ .+ gkm1 .+ 3 .* (hₖ .- hkm1) ./ α
    discr = z .^ 2 .- gₖ .* gkm1; denom = gₖ .+ gkm1 + 2 .* z
    if true in (discr .> 0.0) && true in (abs.(denom) .> eps(Float64))
      #si on peut on utilise l'interpolation cubique
      w = sqrt.(discr)
      dN = -s .* (gₖ .+ z .+ sign.(α) .* w) ./ (denom)
    else #on se rabat sur une étape de sécante
      dN = -gₖ .* s ./ y
    end

    if true in ((tₚ .- t) .* dN .> 0.0) && true in (dN ./ (tₚ .- t) .< γ)
      tplus = t .+ dN
      hplus = obj(h, tplus)
      gplus = grad(h, tplus)
      # verbose && println("N")
    else
      tplus = (t + tₚ) / 2
      hplus = obj(h, tplus)
      gplus = grad(h, tplus)
      # verbose && println("B")
    end

    if true in (t .> tₚ)
      if true in (gplus .< 0.0)
        tₚ = t; tqnp = t; t = tplus
      else
        tqnp = t; t = tplus
      end
    else
      if true in (gplus .> 0.0)
        tₚ = t; tqnp = t; t = tplus
      else
        tqnp = t; t = tplus
      end
    end

    # update values
    hkm1 = hₖ; gkm1 = gₖ; hₖ = hplus; gₖ = gplus
    OK = update_and_stop!(nlpstop, x = t, fx = hₖ, gx = gₖ)
    iter = iter + 1
    verbose &&
        @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                iter, tₚ[1], t[1], dN[1], gₖ[1], gplus[1])
  end

  optimal = nlpstop.meta.optimal
  return optimal, nlpstop
end # function
