export trouve_intervalleA
function trouve_intervalleA(h :: AbstractNLPModel;
                            t₀ :: Float64 = h.meta.x0[1],
                            tₘ :: Float64 = 100.0,
                            direction :: String  =  "Nwt",
                            tol :: Float64 = 1e-10,
                            verbose :: Bool = false,
                            kwargs...)
        # println("on est dans trouve_intervalleA")
        tim1 = t₀
        ti = (tim1 + tₘ) / 2

        i = 0

        h₀ = obj(h, [t₀])[1]; dh₀ = grad(h, [t₀])[1]

        him1 = obj(h, [tim1])[1]; dhim1 = grad(h, [tim1])[1]
        hi = obj(h, [ti])[1]; dhi = grad(h, [ti])[1]

        # println("on est avant les verbose")

        verbose &&
            @printf("iter tim1        dhim1        him1         ti        ")
        verbose &&
            @printf("dhi        hi\n")
        verbose &&
            @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n",
                    i,  tim1, dhim1, him1, ti, dhi, hi)
        # println("avant le while")
        while i < 50
          # print_with_color(:geen, "dans le while i = $i \n")
          him1 = hi
          hi = obj(h, [ti])[1]
          if (hi > h₀) ||((hi  > him1) && (i > 1))
            (topt, iter) = zoom_generic(h, tim1, ti;
                                        direction = direction,
                                        verbose = verbose)
            gopt = grad(h, [topt])[1]
            status = :NotSolved
            (abs(gopt) < 1e-07) && (status = :Optimal)
            (i >= 50) && (status = :Tired)
            tired = i > 50
            optimal = abs(gopt) < 1e-07
            return (topt, obj(h, [topt])[1], norm(gopt, Inf), i, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
          end

          dhim1 = dhi
          dhi = grad(h, [ti])[1]
          if (abs(dhi) <= tol)
            iter = i
            topt = ti
            gopt = grad(h, [topt])[1]
            status = :NotSolved
            (abs(gopt) < 1e-07) && (status = :Optimal)
            (i >= 50) && (status = :Tired)
            tired = i > 50
            optimal = abs(gopt) < 1e-07
            return (topt, obj(h, [topt])[1], norm(gopt, Inf), i, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
          end

          if (dhi >= 0.0)
            (topt, iter) = zoom_generic(h, ti, tim1,
                                        verbose = verbose,
                                        direction = direction)
            gopt = grad(h, [topt])[1]
            status = :NotSolved
            (abs(gopt) < 1e-07) && (status = :Optimal)
            (i >= 50) && (status = :Tired)
            tired = i > 50
            optimal = abs(gopt) < 1e-07
            return (topt, obj(h, [topt])[1], norm(gopt, Inf), i, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)

          end

          tim1 = ti
          ti = (tim1 + tₘ) / 2

          i += 1

          verbose &&
            @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n",
                    i,  tim1, dhim1, him1, ti, dhi, hi)
        end

        return (ti, obj(h, [ti])[1], norm(grad(h, [ti])[1], Inf), i, false, true, :stalled, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
