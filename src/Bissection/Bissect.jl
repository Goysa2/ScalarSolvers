export bissect
function bissect(h       :: AbstractNLPModel,
                      nlpstop :: AbstractStopping;
                      t1 = h.meta.x0, t2 = 100.0,
                      verbose = true)

    (tₐ, tᵦ) = trouve_intervalle(h, t1, 1.0)
    gₐ = grad(h, tₐ)
    gᵦ = grad(h, tᵦ)
    iter = 0

    OK = update_and_start!(nlpstop, x = tₐ, gx = gₐ, g0 = copy(gₐ))

    verbose &&
        @printf(" iter        tₐ        tᵦ         gₐ        gᵦ        \n")
    verbose &&
        @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter, tₐ[1], tᵦ[1], gₐ[1], gᵦ[1])
    while !OK
      tₚ = (tₐ .+ tᵦ) ./ 2
      gₚ = grad(h, tₚ)
      if true in (gₚ .<= 0.0)
        tₐ = tₚ
        gₐ = gₚ
        OK = update_and_stop!(nlpstop, x = tₐ, gx = gₐ)
      else
        tᵦ = tₚ
        gᵦ = gₚ
        OK = update_and_stop!(nlpstop, x = tᵦ, gx = gᵦ)
      end

      iter += 1
      verbose &&
        @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter, tₐ[1], tᵦ[1], gₐ[1], gᵦ[1])
    end

    topt = (tₐ .+ tᵦ) ./ 2.0
    update!(nlpstop.current_state, x = topt, gx = grad(h, topt))
    optimal = OK
    return optimal, nlpstop
end
