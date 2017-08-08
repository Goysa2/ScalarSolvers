export bissect_secA
function bissect_secA(h :: LineModel,
                     t1 :: Float64,
                     t2 :: Float64;
                     tol :: Float64 = 1e-7,
                     maxiter :: Int = 50,
                     verbose :: Bool = false)

      (t0, t1) = trouve_intervalle(h, t1, 2.5)
      γ = 0.8; t = t1; tp = t0; tqnp = t0
      hk = 0.0; hkm1 = 0.0; gkm1 = 0.0; hplus = 0.0
      iter = 0

      hk = obj(h, t); gk = grad(h, t)
      hkm1 = obj(h, tqnp); gkm1 = grad(h, tqnp)

      verbose &&
        @printf(" iter        tp        t         dN         gk          ")
        @printf(" gplus        \n")
      verbose &&
        @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                  iter, tp, t, 0.0, gk, 0.0)

      while ((abs(gk) > tol) && (iter < maxiter)) || (iter == 0)
          s = t - tqnp
          y = gk - gkm1

          Γ = 3 * (gk + gkm1) * s - 6 * (hk - hkm1)
          if y * s + Γ < eps(Float64) * (s^2)
            yt = y
          else
            yt = y + Γ / s
          end
          dN = -gk * s / yt

          if ((tp - t) * dN > 0) && (dN / (tp - t) < γ)
            tplus = t + dN
            hplus = obj(h, tplus)
            gplus = grad(h, tplus)
            verbose && println("N")
          else
            tplus = (t + tp) / 2
            hplus = obj(h, tplus)
            gplus = grad(h, tplus)
            verbose && println("B")
          end

          if t > tp
            if gplus < 0.0
              tp = t; tqnp = t; t = tplus
            else
              tqnp = t; t = tplus
            end
          else
            if gplus > 0.0
              tp = t; tqnp = t; t = tplus
            else
              tqnp = t; t = tplus
            end
          end

          #mise à jour des valeurs
          hkm1 = hk; gkm1 = gk; hk = hplus; gk = gplus
          iter = iter + 1

          verbose &&
            @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n",
                    iter, tp, t, dN, gk, gplus)
        end

      topt = t
      return (topt, iter)
end
