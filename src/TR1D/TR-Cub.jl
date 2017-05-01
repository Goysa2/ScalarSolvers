export TR_Cub
function TR_Cub(h :: AbstractLineFunction,
               t₀ :: Float64,
               tₘ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true,
                eps1 ::Float64=0.2,
                eps2 :: Float64=0.8,
                red :: Float64=0.5,
                aug :: Float64=2.0,
                Δ :: Float64=1.0)

    t = t₀

    iter = 0;
    fₖ = obj(h,t)
    gₖ = grad(h, t)

    secₖ = 1.0
    dN=-gₖ #pour la première itération
    A=0.5
    B=0.0

    q(d) = fₖ + gₖ*d + A*d^2 + B*d^3
    #print_with_color(:yellow," juste avant la boucle while tout va bien \n")

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)
        #print_with_color(:yellow,"on commence la nouvelle itération \n")

        if q(Δ)<q(-Δ)
            d=Δ
        else
            d=-Δ
        end

        #println("avec Δ d=",d)

        cub(t)= fₖ + gₖ*t + A*t^2 + B*t^3
        #dcub(t)= gₖ+ 2*A*t + 3*B*t^2
        #println("on a la dérivée de la cubique")

        #print_with_color(:red,"on a cub et dcub \n")
        #println("fₖ=",fₖ," gₖ=",gₖ," A=",A," B=",B)

        #x₁=(-2*A-sqrt(complex((2*A)^2-4*3*B*gₖ)))/(2*3*B)
        #x₂=(-2*A+sqrt(complex((2*A)^2-4*3*B*gₖ)))/(2*3*B)

        #dR=PolynomialRoots.roots([gₖ,2*A,3*B])
        #dR=PolynomialRoots.roots([-1,-3,1])
        #dR=[x₁ x₂]
        #println("dR=",dR)
        #println("type de dR=",typeof(dR))
        #println("dR[2]=",dR[2])
        #println("dR[1]=",dR[1])
        #println("longeur du tableau de dérivé:")

        #print_with_color(:yellow,"on est capable d'aller calculer dR les racines \n")

        # if isreal(dR[1])
        #   println("on a le premier if dR[1]==Real")
        # end

        #if isreal(dR[2])
        #  println("on a le premier if dR[2]==Real")
        #end

        # if isreal(dR[1]) || isreal(dR[2])
        #   println("on a le premier if dR[1]==Real || dR[2]==Real")
        #   if isreal(dR[1]) & isreal(dR[2])
        #     if cub(dR[1])>cub(dR[2])
        #       dN2=real(dR[2])
        #     else
        #       dN2=real(dR[1])
        #     end
        #   elseif isreal(dR[1]) & !(isreal(dR[2]))
        #     dN2=real(dR[1])
        #   elseif !(isreal(dR[1])) & isreal(dR[2])
        #     dN2=real(dR[2])
        #   end
        # else
        #   dN2 = -gₖ*s/y
        # end
        #println("on a vérifier le isreal")

        #println("gₖ=",gₖ," A=",A," B=",B)

        discr=(2*A)^2-4*3*B*gₖ
        #println("discr=",discr)
        if iter==0 || B==0.0
          dN2=-gₖ/(2*A)
          #println("iter==0 dN2=",dN2)
        else
          if discr>0 #& iter>0
            dR1=(-2*A-sqrt(discr))/(2*3*B)
            dR2=(-2*A+sqrt(discr))/(2*3*B)
            #println("dR1=",dR1," dR2=",dR2)
            if cub(dR1)<cub(dR2)
              dN2=dR1
            else
              dN2=dR2
            end
            #println("avec discre positif et apres comparaison dN2=",dN2)
          elseif discr==0
            dN2=(-2*A)/(2*3*B)
            #println("avec discr==0 dN2=",dN2)
          else
            dN2= -gₖ*s/y
            #println("avec discr<0 dN2=",dN2)
          end
        end

        #print_with_color(:green,"on a dN2 \n")

        if (abs(dN2)<Δ) & (q(d)>q(dN2))
            d=dN2
            #println("on choisit d=dN2")
        end
        #println("d=",d)

        #println("d=",d)

        #print_with_color(:green,"on a d \n")
        #println("d=",d)

        # Numerical reduction computation
        ftestTR = obj(h,t+d)
        gtestTR = grad(h,t+d)

        pred = gₖ*d + A*d^2 + B*d^3

        # truc numérique pour éviter les erreurs d'arrondi lorsque la réduction est faible
        if pred > - 1e-10
            ared = (gₖ + gtestTR)*d/2
        else
            ared = ftestTR-fₖ
        end

        #print_with_color(:yellow,"on a pred et ared \n")

        ratio = ared / pred
        if (ratio < eps1) & (abs(d)==Δ)
            Δ = red*Δ
        else
          #two points memory
          tpred=t
          gkm1=gₖ
          fkm1=fₖ

          t=t+d
          gₖ=gtestTR
          fₖ=ftestTR

          #println("t=",t," tpred=",tpred," gₖ=",gₖ," fₖ=",fₖ," gkm1=",gkm1," fkm1=",fkm1)

          #cubic two points interpolation
          s=t-tpred
          y=gₖ-gkm1

          α=-s
          z= gₖ+gkm1+3*(fₖ-fkm1)/α
          discr=z^2-gₖ*gkm1
          denom=gₖ+gkm1+2*z
          B= 1/3*(gₖ+gkm1+2*z)/(α*α)
          A=-(gₖ+z)/α

          #print_with_color(:green,"on a tout nos paramètres \n")
          #println("A=",A," B=",B)

          if ratio>eps2
            Δ=aug*Δ
          end
        end

        #print_with_color(:green,"on a fini l'itération \n")

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
