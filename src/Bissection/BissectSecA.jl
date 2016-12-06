export bissect_secA
function bissect_secA(h :: AbstractLineFunction,
                    t0 :: Float64,
                    t1 :: Float64;
                    tol :: Float64=1e-7,
                    maxiter :: Int=50,
                    verbose :: Bool=false)

        println("t0=",t0)
        println("t1=",t1)


        if t0==t1
          topt=t0
          iter=0
          return (topt,iter)
        else

          γ=0.8
          t=t1
          println("t=",t)
          tp=t0
          println("tp=",tp)
          tqnp=t0
          hk=0
          hkm1=0
          gkm1=0
          hplus=0
          iter=0

          hk=obj(h,t)
          gk=grad(h,t)
          hkm1=obj(h,tqnp)
          gkm1=grad(h,tqnp)

          verbose && @printf(" iter        tp        t         dN         gk          gplus        \n")
          verbose && @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter,tp,t,0.0,gk, 0.0)

          while ((abs(gk)>tol) & (iter<maxiter)) || (iter==0)
              #println("on entre dans le while: ", iter)
              s = t-tqnp
              #verbose && println("s=",s)
              y = gk-gkm1
              #println("y=",y)

              Γ=3*(gk+gkm1)*s-6*(hk-hkm1)
              verbose && println("gamma=",Γ)
              if y*s+Γ < eps(Float64)*(s^2)
                yt=y
                verbose && println("yt=y")
                verbose && println("yt=",yt)
              else
                yt=y+Γ/s
                verbose && println("yt=y+Γ/s")
                verbose && println("yt=",yt)
              end
              dN=-gk*s/yt
              verbose && println("dN=",dN)

              if ((tp-t)*dN>0) & (dN/(tp-t)<γ)
                tplus = t + dN
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
              iter=iter+1

              verbose && @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter,tp,t,dN,gk,gplus)
            end

          topt=t
          return (topt,iter)
        end


end
