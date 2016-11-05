function ARC_SecA(h :: C1LineFunction,
                t₀ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true)
    #print("on entre dans ")
    # Trust region parameters
    eps1 = 0.1
    eps2 = 0.8
    red = 0.2
    aug = 5
    Δ = 0.5
    t = t₀

    iter = 0;
    hₖ = obj(h,t)
    gₖ = grad(h, t)


    secₖ = 1.0

    q(d) = hₖ + gₖ*d + 0.5*secₖ*d^2

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)

        #step computation
        discr=secₖ^2-4*(gₖ/Δ)
        if discr<0
          discr=secₖ^2+4*(gₖ/Δ)
        end

      #println("on a discr")

        dNp=(-secₖ+sqrt(discr))/(2/Δ) #direction de Newton
        dNp=-2*gₖ/(secₖ+sqrt(discr))

        dNn=(-secₖ-sqrt(discr))/(2/Δ) #direction de Newton
        #println("on a dNn")

        if q(dNp)<q(dNn)
            d=dNp
        else
            d=dNn
        end

        # Numerical reduction computation
        htestTR = obj(h,t+d)
        gtestTR = grad(h,t+d)

        pred = gₖ*d + 0.5*secₖ*d^2

        # truc numérique pour éviter les erreurs d'arrondi lorsque la réduction est faible
        if pred > - 1e-10
            ared = (gₖ + gtestTR)*d/2
        else
            ared = htestTR-hₖ
        end

        ratio = ared / pred
        if (ratio < eps1)
            Δ = red*Δ
        else
            #two point memory
            tpred=t
            gkm1=gₖ
            hkm1=hₖ

            t=t+d
            gₖ=gtestTR
            hₖ=htestTR

            #secant improved two points interpolation
            s=t-tpred
            y=gₖ-gkm1

            Γ=3*(gₖ+gkm1)*s-6*(hₖ-hkm1)
            if y*s+Γ < eps(Float64)*(s^2)
              yt=y
            else
              yt=y+Γ/s
            end

            secₖ=yt/s
            if ratio > eps2
                Δ = aug * Δ
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
