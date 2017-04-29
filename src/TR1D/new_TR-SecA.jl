export new_TR_SecA
function new_TR_SecA(h :: AbstractLineFunction,
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
    fₖ = obj(h,t)
    gₖ = grad(h, t)

    secₖ = 1.0

    q(d) = fₖ + gₖ*d + 0.5*secₖ*d^2

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)

        dS = -gₖ/secₖ; # direction de secante

        d=TR_step_computation(hₖ,secₖ,dS,Δ)

        # Numerical reduction computation
        ftestTR = obj(h,t+d)
        gtestTR = grad(h,t+d)

        pred = gₖ*d + 0.5*secₖ*d^2

        # truc numérique pour éviter les erreurs d'arrondi lorsque la réduction est faible
        if pred > - 1e-10
            ared = (gₖ + gtestTR)*d/2
        else
            ared = ftestTR-fₖ
        end

        ratio = ared / pred
        if (ratio < eps1) & (abs(d)==Δ)
            Δ = red*Δ
        else
            # two points memory
            tpred=t
            gkm1=gₖ
            fkm1=fₖ

            t=t+d
            gₖ=gtestTR
            fₖ=ftestTR

            #secant improved two points interpolation
            s=t-tpred
            y=gₖ-gkm1

            Γ=3*(gₖ+gkm1)*s-6*(fₖ-fkm1)
            if (y*s+Γ)<eps(Float64)*(s^2) #correction trop petite
              yt=y
            else
              yt=y+Γ/s
            end

            secₖ=yt/s

            if ratio>eps2
              Δ=aug *Δ
            end
          end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
