export Trouve_IntervalleA
function Trouve_IntervalleA(h :: AbstractNLPModel,
                            nlpstop :: AbstractStopping;
                            t₀ = h.meta.x0, tₘ = 100.0,
                            direction =  :Nwt, verbose = false, kwargs...)

  tim1 = t₀
  ti = (tim1 .+ tₘ) ./ 2

  i = 0

  h₀ = obj(h, t₀); dh₀ = grad(h, t₀)

  him1 = obj(h, tim1); dhim1 = grad(h, tim1)
  hi = obj(h, ti); dhi = grad(h, ti)

  OK = update_and_start!(nlpstop, x = ti, fx = hi, gx = dhi, g0 = copy(dh₀))

  verbose &&
    @printf("iter tim1        dhim1        him1         ti        ")
  verbose &&
    @printf("dhi        hi\n")
  verbose &&
    @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i,  tim1[1], dhim1[1], him1[1], ti[1], dhi[1], hi[1])

  while !OK
    him1 = hi
    hi = obj(h, ti)
    update!(nlpstop.current_state, fx = hi)
    if true in (hi .> h₀) || true in ((hi  .> him1) && (i > 1))
      optimal, nlpstop  = zoom_generic(h, nlpstop, tim1, ti;
                                       direction = direction, verbose = verbose)
      return optimal, nlpstop
    end

    dhim1 = dhi; dhi = grad(h, ti)
    update!(nlpstop.current_state, gx = dhi)

    if (dhi .>= 0.0)
      optimal, nlpstop = zoom_generic(h, nlpstop, ti, tim1, verbose = verbose, direction = direction)
      return optimal, nlpstop
    end

    tim1 = ti
    ti = (tim1 .+ tₘ) ./ 2

    OK = update_and_stop!(nlpstop, x = ti)

    i += 1

    verbose &&
      @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i,  tim1[1], dhim1[1], him1[1], ti[1], dhi[1], hi[1])
  end

  return OK, nlpstop
end
