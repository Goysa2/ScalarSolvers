export bissect
function bissect(h :: AbstractNLPModel;
                t1 :: Float64 = -10.0,
                t2 :: Float64 = 100.0,
                tol :: Float64 = 1e-7,
                maxiter :: Int = 50,
                verbose :: Bool = true)

    (tₐ, tᵦ) = trouve_intervalle(h, t1, 1.0)
    gₐ = grad(h, [tₐ])[1]
    gᵦ = grad(h, [tᵦ])[1]
    iter = 0

    verbose &&
        @printf(" iter        tₐ        tᵦ         gₐ        gᵦ        \n")
    verbose &&
        @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter, tₐ, tᵦ, gₐ, gᵦ)
    while ((max(abs(gᵦ), abs(gₐ))) > tol) && (iter < maxiter) || iter == 0
      tₚ = (tₐ + tᵦ) / 2
      gₚ = grad(h, [tₚ])[1]
      if gₚ <= 0.0
        tₐ = tₚ
        gₐ = gₚ
      else
        tᵦ = tₚ
        gᵦ = gₚ
      end

      iter += 1
      verbose &&
        @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter, tₐ, tᵦ, gₐ, gᵦ)
    end

    topt = (tₐ + tᵦ) / 2.0
    gopt = grad(h, [topt])[1]
    if maxiter <= iter
        tired = true
    else
        tired = false
    end
    if (abs(gopt) > tol)
        optimal = false
    else
        optimal = true
    end
    status = :tmp
    return (topt, obj(h, [topt])[1], norm(gopt, Inf), iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
end
