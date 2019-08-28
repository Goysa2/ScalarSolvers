export bissect_secA
function bissect_secA(h :: AbstractNLPModel,
                      nlpstop :: AbstractStopping;
                      t1 = h.meta.x0, t2 = 100.0,
                      verbose :: Bool = false)


    (t0, t1) = trouve_intervalle(h, t1, 2.5)
    γ = 0.8; t = t1; tp = t0; tqnp = t0
    hk = [0.0]; hkm1 = [0.0]; gkm1 = [0.0]; hplus = [0.0];
    iter = 0

    hk   = obj(h, t);    gk   = grad(h, t)
    hkm1 = obj(h, tqnp); gkm1 = grad(h, tqnp)

    OK = update_and_start!(nlpstop, x = t, gx = gk, g0 = copy(gk))

    verbose &&
      @printf(" iter        tp        t         dN         gk          ")
    verbose &&
      @printf(" gplus        \n")
    verbose &&
      @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                iter, tp[1], t[1], 0.0, gk[1], 0.0)

    while !OK
        s = t .- tqnp
        y = gk .- gkm1

        Γ = [3.0] .* (gk .+ gkm1) .* s .- [6.0] .* (hk .- hkm1)
        if y .* s .+ Γ < eps.(Float64) .* (s .^ 2)
          yt = y
        else
          yt = y .+ Γ ./ s
        end
        dN = vec(-gk .* s ./ yt)

        if  true in ((tp - t) .* dN .> 0.0) &&  true in (dN ./ (tp - t) .< γ)
          tplus = t .+ dN
          hplus = obj(h, tplus)
          gplus = grad(h, tplus)
        else
          tplus = (t .+ tp) ./ 2
          hplus = obj(h, tplus)
          gplus = grad(h, tplus)
        end

        if (true in (t .> tp))
          if (true in (gplus .< 0.0))
            tp = t; tqnp = t; t = tplus
          else
            tqnp = t; t = tplus
          end
        else
          if (true in (gplus .> 0.0))
            tp = t; tqnp = t; t = tplus
          else
            tqnp = t; t = tplus
          end
        end

        # update values
        hkm1 = hk; gkm1 = gk; hk = hplus; gk = gplus
        OK = update_and_stop!(nlpstop, x = tp, gx = gk)
        iter = iter + 1

        verbose &&
          @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                  iter, tp[1], t[1], dN[1], gk[1], gplus[1])
      end

      optimal = OK
      return optimal, nlpstop
end
