export bissect_secA
function bissect_secA(h :: AbstractNLPModel;
                     t1 :: Float64 = h.meta.x0[1],
                     t2 :: Float64 = 100.0,
                     tol :: Float64 = 1e-7,
                     maxiter :: Int = 50,
                     verbose :: Bool = false)

    (length(h.meta.x0) > 1) && warn("Not a 1-D problem ")
      (t0, t1) = trouve_intervalle(h, t1, 2.5)
      γ = 0.8; t = t1; tp = t0; tqnp = t0
      hk = 0.0; hkm1 = 0.0; gkm1 = 0.0; hplus = 0.0
      iter = 0

      hk = obj(h, [t])[1]; gk = grad(h, [t])[1]
      hkm1 = obj(h, [tqnp])[1]; gkm1 = grad(h, [tqnp])[1]

      verbose &&
        @printf(" iter        tp        t         dN         gk          ")
      verbose &&
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
            hplus = obj(h, [tplus])[1]
            gplus = grad(h, [tplus])[1]
            verbose && println("N")
          else
            tplus = (t + tp) / 2
            hplus = obj(h, [tplus])[1]
            gplus = grad(h, [tplus])[1]
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
      status = :NotSolved
      (abs(gk) < tol) && (status = :Optimal)
      (iter >= maxiter) && (status = :Tired)
      tired = iter > maxiter
      optimal = abs(gk) < tol
      return (topt, obj(h, [topt])[1], norm(gk, Inf), iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
