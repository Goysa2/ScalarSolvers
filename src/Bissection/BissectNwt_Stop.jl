export bissect_nwt_stop
function bissect_nwt_stop(h       :: AbstractNLPModel,
                          nlpstop :: AbstractStopping;
                          t1 = h.meta.x0, t2 = 100.0,
                          tol = 1e-7, maxiter = 50, verbose = false)

  (tₐ, tᵦ) = trouve_intervalle_stop(h, t1, 2.5)
  γ = 0.8; t = tᵦ; tₚ = tₐ;
  hₖ = [0.0]; hplus = [0.0]
  iter = 0

  gₖ = grad(h, t)

  OK = update_and_start!(nlpstop, x = t, gx = gₖ, g0 = copy(gₖ))

  verbose &&
    @printf(" iter        tₚ        t         dN         gₖ          ")
  verbose &&
    @printf("gplus        \n")
  verbose &&
    @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                iter, tₚ[1], t[1], 0.0, gₖ[1], 0.0)
  # printstyled("on est juste avant le while \n", color = :bold)
  while !OK

      kₖ = hess(h, t)
      dN = vec(-gₖ ./ kₖ) #direction de Newton
      if (true in ((tₚ .- t) .* dN .> 0.0)) && (true in (dN ./ (tₚ - t) .< γ))
        tplus = t .+ dN
        #hplus = obj(h, tplus)
        gplus = grad(h, tplus)
        # verbose && println("N")
      else
        tplus = (t .+ tₚ) ./ 2
        #hplus = obj(h, tplus)
        gplus = grad(h, tplus)
        # verbose && println("B")
      end

      if true in (t .> tₚ)
        if true in (gplus .< 0.0)
          tₚ = t;  t = tplus
        else
           t = tplus
        end
      else
        if true in (gplus .> 0.0)
          tₚ = t;  t = tplus
        else
          t = tplus
        end
      end

      #mise à jour des valeurs
      hₖ = hplus; gₖ = gplus
      OK = update_and_stop!(nlpstop, x = tₚ, gx = gₖ)
      iter = iter + 1
      verbose &&
        @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                iter, tₚ[1], t[1], dN[1], gₖ[1], gplus[1])
    end

  optimal = OK
  return optimal, nlpstop
end
