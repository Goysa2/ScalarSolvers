export bissect
function bissect(h :: C2LineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true)

        g=grad(h,t₀)
        k=hess(h,t₀)
        inc=abs(g/k)

        (tₐ,tᵦ,iter)=trouve_intervalle(h,t₀,inc)

        #tₐ=t₀
        #tᵦ=t₁
        hₐ=obj(h,tₐ)
        gₐ=grad(h,tₐ)
        hᵦ=obj(h,tᵦ)
        gᵦ=grad(h,tᵦ)
        iter=0

        #println("on défini les paramètres de la bissection")

        verbose && @printf(" iter        tₐ        tᵦ         gₐ        gᵦ        \n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter,tₐ,tᵦ,gₐ,gᵦ)
        #println("on est juste avant la boucle while")
        while ((max(abs(gᵦ),abs(gₐ)))>tol) & (iter<maxiter) || iter==0
          #println("on est entré dans la boucle while")
          tₚ=(tₐ+tᵦ)/2
          #mise à jour des valeur
          gₚ=grad(h,tₚ)
          if gₚ<=0
            tₐ=tₚ
            gₐ=gₚ
          else
            tᵦ=tₚ
            gᵦ=gₚ
          end

          iter=iter+1
          verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter,tₐ,tᵦ,gₐ,gᵦ)
        end

        topt=(tₐ+tᵦ)/2

        return (topt, iter)

end
