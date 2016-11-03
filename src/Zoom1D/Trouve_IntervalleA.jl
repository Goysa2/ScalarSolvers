function trouve_intervalleA(h :: C2LineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                ϵ :: Float64=1e-10,
                verbose :: Bool=true)

        tim1=t₀
        ti=(tim1+tₘ)/2

        i=0

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        him1=obj(h,tim1)
        dhim1=grad(h,tim1)
        hi=obj(h,ti)
        dhi=grad(h,ti)

        #println("on les paramètre de trouve_intervalleA")

        verbose && @printf("iter tim1        dhim1        him1         ti        dhi        hi\n")
        verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)

        if (abs(dhi)<ϵ)
          topt=ti
          iter=i
          return (topt,iter)
        end

        while  i<50
          #println("On est dans le while de trouve_intervalleA")
          if ((hi>h₀+0.01*dh₀) ||(hi>him1)) & (i>1)
            #println("1er if")
            (topt,iter)=zoom(h,tim1,ti)
            if obj(h,topt)<obj(h,t₀)
              println("h(0)<h(t*)")
            end
            return (topt,iter)
          end

          if (abs(dhi)<=ϵ)
            #println("2e if")
            iter=i
            topt=ti
            # if obj(h,topt)<obj(h,t₀)
            #   println("h(0)<h(t*)")
            # end
            return (topt,iter)
          end

          if (dhi>=0)
            #println("3e if")
            (topt,iter)=zoom(h,ti,tim1)
            # if obj(h,topt)<obj(h,t₀)
            #   println("h(t*)<h(t₀)")
            # end
            return (topt,iter)
          end

          tim1=ti
          ti=(tim1+tₘ)/2

          him1=obj(h,tim1)
          dhim1=grad(h,tim1)
          hi=obj(h,ti)
          dhi=grad(h,ti)



          i=i+1
          #println("on continue d'itérer")

          verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)
        end

        # if obj(h,topt)<obj(h,t₀)
        #   println("h(0)<h(t*)")
        # end
end
