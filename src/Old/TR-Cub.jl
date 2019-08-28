export TR_Cub
function TR_Cub(h :: AbstractNLPModel;
               t₀ :: Float64 = h.meta.x0[1],
               tol :: Float64 = 1e-7,
               maxiter :: Int = 50,
               verbose :: Bool = true,
               eps1 ::Float64 = 0.2,
               eps2 :: Float64 = 0.8,
               red :: Float64 = 0.5,
               aug :: Float64 = 2.0,
               Δ :: Float64 = 1.0)

    (length(h.meta.x0) > 1) && warn("Not a 1-D problem ")
    t = t₀; iter = 0;
    fₖ = obj(h, [t])[1]
    gₖ = grad(h, [t])[1]

    secₖ = 1.0
    dN = -gₖ #pour la première itération
    A = 0.5; B = 0.0

    q(d) = fₖ + gₖ * d + A * d^2 + B * d^3

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ) > tol) & (iter < maxiter)) | (iter == 0)

        if q(Δ) < q(-Δ)
            d = Δ
        else
            d = -Δ
        end

        cub(t) = fₖ + gₖ * t + A * t^2 + B * t^3
        #dcub(t)=gₖ+2*A*t+3*B*t^2
        dR = roots(Poly([gₖ, 2 * A, 3 * B]))

        if isreal(dR)
          dR = real(dR)
          dN2 = dR[1]
          if length(dR) > 1
            if (cub(dR[2]) < cub(dR[1]))
              dN2 = dR[2]
            end
          end
        else
          dN2 = -gₖ * s / y
        end

        if (abs(dN2) < Δ) & (q(d) > q(dN2))
            d = dN2
        end

        if d == 0.0
            return (t, fₖ, norm(gₖ, Inf), iter + 1, false, false, :FailedDirectionComputation, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)
            #L'OUTIL DE CALCUL DES RACINES DU PACKAGE Polynomials ARRONDIE À 0 CERTAINES RACINES D'OÙ LA NOUVELLE CONDITION
        end

        ftestTR = obj(h, [t + d])[1]
        gtestTR = grad(h, [t + d])[1]

        pred = gₖ * d + A * d^2 + B * d^3

        if pred > - 1e-10
            ared = (gₖ + gtestTR) * d / 2
        else
            ared = ftestTR - fₖ
        end

        ratio = ared / pred
        if (ratio < eps1) & (abs(d) == Δ)
            Δ = red * Δ
        else
          #two points memory
          tpred = t; gkm1 = gₖ; fkm1 = fₖ
          t = t + d; gₖ = gtestTR; fₖ = ftestTR

          #cubic two points interpolation
          s = t - tpred; y = gₖ - gkm1

          α = -s
          z = gₖ + gkm1 + 3 * (fₖ - fkm1) / α
          discr = z^2 - gₖ * gkm1
          denom = gₖ + gkm1 + 2 * z
          B = 1 / 3 * (gₖ + gkm1 +2 * z) / (α * α)
          A = -(gₖ + z) / α

          if ratio > eps2
            Δ = aug * Δ
          end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    status = :NotSolved
    (abs(gₖ) < tol) && (status = :Optimal)
    (iter >= maxiter) && (status = :Tired)
    tired = iter > maxiter
    optimal = abs(gₖ) < tol
    return (t, fₖ, norm(gₖ, Inf), iter, optimal, tired, status, h.counters.neval_obj, h.counters.neval_grad, h.counters.neval_hess)

end
