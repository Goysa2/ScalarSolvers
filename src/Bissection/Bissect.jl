export bissect
function bissect(h :: AbstractLineFunction2,
                tₐ :: Float64,
                tᵦ :: Float64;
                tol :: Float64=1e-7,
                maxiter :: Int=50,
                verbose :: Bool=true)

      if tₐ==tᵦ
        topt=tᵦ
        iter=0
        return (topt,iter)
      else
        gₐ=grad(h,tₐ)
        gᵦ=grad(h,tᵦ)
        iter=0

        verbose && @printf(" iter        tₐ        tᵦ         gₐ        gᵦ        \n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e\n", iter,tₐ,tᵦ,gₐ,gᵦ)
        while ((max(abs(gᵦ),abs(gₐ)))>tol) & (iter<maxiter) || iter==0
          tₚ=(tₐ+tᵦ)/2
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
end
