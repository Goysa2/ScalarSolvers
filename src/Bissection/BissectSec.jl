export bissect_sec
function bissect_sec(h :: AbstractLineFunction2,
                     tₐ :: Float64,
                     tᵦ :: Float64;
                     tol :: Float64=1e-7,
                     maxiter :: Int=50,
                     verbose :: Bool=false)

          if tₐ==tᵦ
           topt=tᵦ
           iter=0
           return (topt,iter)
          else

              γ=0.8
              t=tᵦ
              tₚ=tₐ
              tqnp=tₐ
              hₖ=0
              hkm1=0
              gkm1=0
              hplus=0
              iter=0

              gₖ=grad(h,t)
              gkm1=grad(h, tqnp)

              verbose && @printf(" iter        tₚ        t         dN         gₖ          gplus        \n")
              verbose && @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter,tₚ,t,0.0,gₖ, 0.0)

              while ((abs(gₖ)>tol) & (iter<maxiter)) || (iter==0)
                #println("on entre dans le while: ", iter)
                s = t-tqnp
                y = gₖ-gkm1

                dN=-gₖ*s/y #direction de Newton

                if ((tₚ-t)*dN>0) & (dN/(tₚ-t)<γ)
                  tplus = t + dN
                  #hplus = obj(h, tplus)
                  gplus = grad(h,tplus)
                  verbose && print_with_color(:green,"N")
                else
                  tplus = (t+tₚ)/2
                  #hplus = obj(h, tplus)
                  gplus = grad(h,tplus)
                  verbose && print_with_color(:green,"B")
                end

                if t>tₚ
                  if gplus<0
                    tₚ=t
                    tqnp=t
                    t=tplus
                  else
                    tqnp=t
                    t=tplus
                  end
                else
                  if gplus>0
                    tₚ=t
                    tqnp=t
                    t=tplus
                  else
                    tqnp=t
                    t=tplus
                  end
                end

                #mise à jour des valeurs
                hkm1=hₖ
                gkm1=gₖ
                hₖ=hplus
                gₖ=gplus
                iter=iter+1
                verbose && @printf(" %7.2e %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter,tₚ,t,dN,gₖ,gplus)
              end

                topt=t
                return (topt,iter)

          end
end
