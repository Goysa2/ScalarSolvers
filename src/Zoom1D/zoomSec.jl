export zoom_Sec
function zoom_Sec(h :: AbstractLineFunction2,
                  t₀ :: Float64,
                  t₁ :: Float64;
                  c₁ :: Float64=0.01,
                  tol :: Float64=1e-7,
                  ϵ :: Float64=1e-10,
                  maxiter :: Int=100,
                  verbose :: Bool=false)
                  
        if obj(h,t₀)<obj(h,t₁)
          tl=t₀
          th=t₁
        else
          tl=t₁
          th=t₀
        end

        hl=obj(h,tl)
        dhl=grad(h,tl)
        hh=obj(h,th)
        dhh=grad(h,th)

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        γ=0.8
        t=t₁
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

        verbose && println("on a les paramètre de zoom")

        verbose && @printf(" iter        tₗ        tₕ         t        hₗ        hₕ         hk         gk\n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,t,hl,hh,hk,gk)


        while ((i<maxiter) & (abs(gk)>tol)) || i==0
          if i==0
            verbose && println("on est dans le while de zoom")
          end
          #hk=obj(h,t)
          if (hk>h₀+c₁*(t-t₀)*dh₀) || (hl<=hk)
            if (hk>h₀+c₁*t*dh₀)
              verbose && println("(hk>h₀+c₁*ti*dh₀)")
            else
              verbose && println("(hl<=hk)")
            end
            th=t
            hh=hk
            dhh=gk
          else
            verbose && verbose && println("else")
            #gk=grad(h,t)
            if (abs(gk)<ϵ)
              verbose && print(abs(gk),"<",-0.99*dh₀)
              verbose && println("1er if dans le else zoom")
              topt=t
              iter=i
              return (topt,iter)
            elseif gk*(th-tl)>=0
              verbose && println("elseif dans le else de zoom")
              th=tl
              hh=hl
              dhh=dhl
            end
            tl=t
            hl=hk
            dhl=gk
          end

          s = t-tqnp
          y = gk-gkm1

          dS=-gk*s/y

          if ((tp-t)*dS>0) & (dS/(tp-t)<γ)
            tplus = t + dS
            hplus = obj(h, tplus)
            gplus = grad(h,tplus)
            verbose && println("N")
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

          #mise à jour des valeurs
          hkm1=hk
          gkm1=gk
          hk=hplus
          gk=gplus
          i=i+1
          #println("on continue d'itérer")
          verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,t,hl,hh,hk,gk)

        end
        topt=t
        iter=i

        return (topt,iter)
end
