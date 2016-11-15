export trouve_intervalleA
function trouve_intervalleA(h :: AbstractLineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                ϵ :: Float64=1e-10,
                verbose :: Bool=false)

        tim1=t₀
        ti=(tim1+tₘ)/2

        i=0

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        him1=obj(h,tim1)
        dhim1=grad(h,tim1)
        hi=obj(h,ti)
        dhi=grad(h,ti)

        verbose && @printf("iter tim1        dhim1        him1         ti        dhi        hi\n")
        verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)

        if (abs(dhi)<ϵ)
          topt=ti
          iter=i
          return (topt,iter)
        end

        while  i<50
          if ((hi>h₀+0.01*dh₀) ||(hi>him1)) & (i>1)
            (topt,iter)=zoom(h,tim1,ti)
            return (topt,iter)
          end

          if (abs(dhi)<=ϵ)
            iter=i
            topt=ti
            return (topt,iter)
          end

          if (dhi>=0)
            (topt,iter)=zoom(h,ti,tim1)
            return (topt,iter)
          end

          tim1=ti
          ti=(tim1+tₘ)/2

          him1=hi
          dhim1=dhi
          hi=obj(h,ti)
          dhi=grad(h,ti)

          i=i+1

          verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)
        end

end
