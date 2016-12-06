export trouve_intervalleACub
function trouve_intervalleACub(h :: AbstractLineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                ϵ :: Float64=1e-10,
                verbose :: Bool=false)
        print_with_color(:magenta,"trouve_intervalleACub")

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        γ=0.8
        t=tₘ
        tp=t₀
        tqnp=t₀
        hk=0
        hkm1=0
        gkm1=0
        hp=0
        i=0

        hk=obj(h,t)
        gk=grad(h,t)
        hkm1=obj(h,tqnp)
        gkm1=grad(h,tqnp)

        verbose && @printf("iter tp        gkm1        hkm1         t        gk        hk\n")
        verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tp,gkm1,hkm1,t,gk,hk)

        while  i<50
          if (hk>h₀+0.01*(t-t₀)*dh₀) ||((hk>hkm1) & (i>1))
            (topt,iter)=zoom_Cub(h,tp,t,verbose=false)
            return (topt,iter)
          end

          if (abs(gk)<=ϵ)
            iter=i
            topt=t
            return (topt,iter)
          end

          if (gk>=0)
            (topt,iter)=zoom_Cub(h,t,tp,verbose=false)
            return (topt,iter)
          end
          s = t-tqnp
          y = gk-gkm1

          α=-s
          z=gk+gkm1+3*(hk-hkm1)/α
          discr=z^2-gk*gkm1
          denom=gk+gkm1+2*z
          if (discr>0) & (abs(denom)>eps(Float64))
            #si on peut on utilise l'interpolation cubique
            w=sqrt(discr)
            dC=-s*(gk+z+sign(α)*w)/(denom)
          else #on se rabat sur une étape de sécante
            dC=-gk*s/y
          end

          if ((tp-t)*dC>0) & (dC/(tp-t)<γ)
            tplus = t + dC
            hplus = obj(h, tplus)
            gplus = grad(h,tplus)
            verbose && println("C")
          else
            tplus = (t+tp)/2
            hplus = obj(h, tplus)
            gplus = grad(h,tplus)
            verbose && println("B")
          end

          if t>tp
            if gplus<0
              tp=t
              tqnp=t
              t=tplus
            else
              tqnp=t
              t=tplus
            end
          else
            if gplus>0
              tp=t
              tqnp=t
              t=tplus
            else
              tqnp=t
              t=tplus
            end
          end

          hkm1=hk
          gkm1=gk
          hk=hplus
          gk=gplus
          i=i+1

          verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,gkm1,hkm1,ti,gk,hk)
        end

end
