export new_TR_generic
function new_TR_generic(h :: AbstractLineFunction2,
                        t₀ :: Float64,
                        tₘ :: Float64;
                        tol :: Float64=1e-7,
                        maxiter :: Int=50,
                        verbose :: Bool=true,
                        eps1 :: Float64=0.2,
                        eps2 :: Float64=0.8,
                        red:: Float64=0.5,
                        aug :: Float64=2.0,
                        Δ :: Float64 = 1.0,
                        direction :: String="Nwt")

    t = t₀

    iter = 0;
    fₖ = obj(h,t)
    gₖ = grad(h, t)

    if direction=="Sec" || direction=="SecA"
      secₖ = 1.0
    else
      hₖ=hess(h,t)
    end

    #q(d) = fₖ + gₖ*d + 0.5*hₖ*d^2 #|| q(d) = fₖ + gₖ*d + 0.5*secₖ*d^2

    verbose && @printf(" iter  t         gₖ          Δ        pred         ared\n")
    verbose && @printf(" %4d %7.2e  %7.2e  %7.2e \n", iter,t,gₖ,Δ)

    while ((abs(gₖ)>tol) & (iter < maxiter)) | (iter == 0)

        if direction=="Sec" || direction=="SecA"
          dS = -gₖ/secₖ; # direction de secante
          d=TR_step_computation(secₖ,gₖ,dS,Δ)
        elseif direction=="Nwt"
          dN = -gₖ/hₖ; # direction de Newton
          d=TR_step_computation(hₖ,gₖ,dN,Δ)
          #println("dN=",dN)
          #println("d=",d)
        end

        # Numerical reduction computation
        ftestTR = obj(h,t+d)
        gtestTR = grad(h,t+d)

        if direction=="Nwt"
          (pred,ared,ratio)=pred_ared_computation(gₖ,fₖ,hₖ,d,ftestTR,gtestTR)
        elseif direction=="Sec" || direction=="SecA"
          (pred,ared,ratio)=pred_ared_computation(gₖ,fₖ,secₖ,d,ftestTR,gtestTR)
        end

        #println("(pred,ared,ratio)=(",pred,",",ared,",",ratio,")")

        if direction=="Nwt"
          cond_red=(ratio < eps1)
        else
          cond_red=(ratio < eps1) & (abs(d)==Δ)
        end

        if cond_red
            Δ = red*Δ
        else
            #println("on est dans le else")
            if direction=="Sec"
              (t,fₖ,gₖ,s,y,secₖ)  =Sec_computation(t,gₖ,d,ftestTR,gtestTR)
            elseif direction=="SecA"
              (t,fₖ,gₖ,secₖ)  =SecA_computation(t,gₖ,fₖ,d,ftestTR,gtestTR)
            elseif direction=="Nwt"
              (t,gₖ,fₖ,hₖ)=Nwt_computation(t,d,gtestTR,ftestTR,h)
              #println("après Nwt_computation (t,gₖ,fₖ,hₖ)=(",t,",",gₖ,",",fₖ,",",hₖ,")")
            end

            if ratio > eps2
                Δ = aug * Δ
            end
        end

        iter += 1
        verbose && @printf(" %4d %7.2e  %7.2e  %7.2e %7.2e %7.2e\n", iter,t,gₖ,Δ,pred,ared)
    end

    return (t, iter)
end
