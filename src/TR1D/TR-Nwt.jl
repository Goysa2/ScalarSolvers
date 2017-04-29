export TR_Nwt
function TR_Nwt(hh :: AbstractLineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true,
                eps1 :: Float64=0.2,
                eps2 :: Float64=0.8,
                red :: Float64=0.5,
                aug :: Float64=2.0,
                Δ :: Float64=1.0)


    t = t₀

    iter = 0;
    fₖ = obj(hh,t)
    gₖ = grad(hh, t)
    hₖ = hess(hh,t)

    q(d) = fₖ + gₖ*d + 0.5*hₖ*d^2

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) || (iter == 0)

        dN = -gₖ/hₖ; # direction de Newton

        if q(Δ)<q(-Δ)
            d=Δ
        else
            d=-Δ
        end

        if (abs(dN)<Δ) & (q(d)>q(dN))
            d=dN
        end

        # Numerical reduction computation
        ftestTR = obj(hh,t+d)
        gtestTR = grad(hh,t+d)

        pred = gₖ*d + 0.5*hₖ*d^2

        # truc numérique pour éviter les erreurs d'arrondi lorsque la réduction est faible
        if pred > - 1e-10
            ared = (gₖ + gtestTR)*d/2
        else
            ared = ftestTR-fₖ
        end

        ratio = ared / pred
        if (ratio < eps1)
            Δ = red*Δ
        else
            t=t+d

            gₖ=gtestTR
            fₖ=ftestTR
            hₖ=hess(hh,t)

            if ratio > eps2
                Δ = aug * Δ
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
