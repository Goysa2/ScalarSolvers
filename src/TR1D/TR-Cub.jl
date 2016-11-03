function TR_Cub(h :: C2LineFunction,
               t₀ :: Float64,
               tₘ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true)
    #print("on entre dans ")
    # Trust region parameters
    eps1 = 0.2
    eps2 = 0.8
    red = 0.5
    aug = 2
    Δ = 1.0
    t = t₀

    iter = 0;
    hₖ = obj(h,t)
    gₖ = grad(h, t)

    secₖ = 1.0
    dN=-gₖ #pour la première itération
    A=0.5
    B=0.0

    q(d) = hₖ + gₖ*d + A*d^2 + B*d^3

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)

        if q(Δ)<q(-Δ)
            d=Δ
        else
            d=-Δ
        end

        cub(t)= hₖ + gₖ*t + A*t^2 + B*t^3
        dcub(t)= gₖ+ 2*A*t + 3*B*t^2
        #println("on a la dérivée de la cubique")


        dR=roots(dcub)
        #println("longeur du tableau de dérivé:")


        if isreal(dR)
          #println("on a le premier if")
          dR=real(dR)
          dN2=dR[1]
          #println("on a le premier dN2")
          if (length(dR)>1)
            if (cub(dR[1])>cub(dR[2]))
            dN2=dR[2]
            #println("on est dans le deuxième if et l'autre dN2")
            end
          end
        else
          dN2 = -gₖ*s/y
        end
        #println("on a vérifier le isreal")

        if (abs(dN2)<Δ) & (q(d)>q(dN2))
            d=dN2
        end

        # Numerical reduction computation
        htestTR = obj(h,t+d)
        gtestTR = grad(h,t+d)

        pred = gₖ*d + A*d^2 + B*d^3

        # truc numérique pour éviter les erreurs d'arrondi lorsque la réduction est faible
        if pred > - 1e-10
            ared = (gₖ + gtestTR)*d/2
        else
            ared = htestTR-hₖ
        end

        ratio = ared / pred
        if (ratio < eps1) & (abs(d)==Δ)
            Δ = red*Δ
        else
          #two points memory
          tpred=t
          gkm1=gₖ
          hkm1=hₖ

          t=t+d
          gₖ=gtestTR
          hₖ=htestTR

          #cubic two points interpolation
          s=t-tpred
          y=gₖ-gkm1

          α=-s
          z= gₖ+gkm1+3*(hₖ-hkm1)/α
          discr=z^2-gₖ*gkm1
          denom=gₖ+gkm1+2*z
          B= 1/3*(gₖ+gkm1+2*z)/(α*α)
          A=-(gₖ+z)/α

          if ratio>eps2
            Δ=aug*Δ
          end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
