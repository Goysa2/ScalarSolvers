export trouve_intervalleASec
function trouve_intervalleASec(h :: C2LineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                ϵ :: Float64=1e-10,
                verbose :: Bool=false)

        ϵ₂=1e-5

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        t=tₘ
        ti=t
        tim1=t₀
        tqnp=t₀
        him1=0
        dhim1=0
        hi=0
        dhi=0
        i=0
        #
        # s=t-tqnp
        # y=dhi-dhim1
        #
        # println("s=",s," y=",y)
        #
        # dS=-dhi*s/y
        #
        # tiB=(tim1+t)/2
        # tiS=t+dS
        #
        # println("tiB=",tiB," tiS=",tiS)
        # i=0
        #
        # if (tiS<t₀) || (tiS>tₘ) || (tiS==NaN)
        #   println("on choisit tiB")
        #   ti=tiB
        #   hi=obj(h,tiB)
        #   dhi=grad(h,tiB)
        #   print_with_color(:green,"B")
        # else
        #   println("on choisit tiS")
        #   ti=tiS
        #   hi=obj(h,tiS)
        #   dhi=grad(h,tiS)
        #   print_with_color(:green,"S")
        # end
        #
        # if t>ti
        #   if dhi<0
        #     tim1=t
        #     tqnp=t
        #     t=ti
        #   else
        #     tqnp=t
        #     t=ti
        #   end
        # else
        #   if dhi>0
        #     tim1=t
        #     tqnp=t
        #     t=ti
        #   else
        #     tqnp=t
        #     t=ti
        #   end
        # end

        #println("on a les paramètres de trouve intervalle")

        verbose && @printf("iter tim1      dhim1      him1       ti      dhi      hi\n")
        verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)

        # if (abs(dhi)<ϵ)
        #   topt=t
        #   iter=i
        #   return (topt,iter)
        # end

        while i==0 || ((abs(dhi)>ϵ) & (i<50))
          # if i==0
          #  println("on est dans le while de trouve intervalle")
          # end

          if ((hi>h₀+0.01*t*dh₀) ||(hi>him1)) & (i>1)
            #println("1er if de trouve_intervalleASec")
            (topt,iter)=zoom_Sec(h,tim1,ti)
            # if obj(h,topt)<obj(h,t₀)
            #  println("h(0)<h(t*)")
            # end
            return (topt,iter)
          end

          if (i>0) & (abs(dhi)<=ϵ)
            #println("2e if de trouve_intervalleASec ")
            iter=i
            topt=ti
            # if obj(h,topt)<obj(h,t₀)
            #   println("h(0)<h(t*)")
            # end
            return (topt,iter)
          end

          #if (dhi*dhim1<=0) #pas la vrai condition...
          if (dhi)>=0
            #println("3e if de trouve_intervalleASec")
            (topt,iter)=zoom_Sec(h,ti,tim1)
            # if obj(h,topt)<obj(h,t₀)
            #   println("h(t*)<h(t₀)")
            # end
            return (topt,iter)
          end

          s=t-tqnp
          y=dhi-dhim1

          #println("s=",s," y=",y)
          dS=-dhi*s/y

          tiB=(ti+t)/2
          tiS=t+dS

          him1=hi
          dhim1=dhi

          #println("tiB=",tiB," tiS=",tiS)

          if (tiS<t₀) || (tiS>tₘ)
            #println("on choisit tiB")
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            verbose && print_with_color(:green,"B")
          else
            #println("on choisit tiS")
            ti=tiS
            hi=obj(h,tiS)
            dhi=grad(h,tiS)
            verbose && print_with_color(:green,"S")
          end

          if t>ti
            if dhi<0
              tim1=t
              tqnp=t
              t=ti
            else
              tqnp=t
              t=ti
            end
          else
            if dhi>0
              tim1=t
              tqnp=t
              t=ti
            else
              tqnp=t
              t=ti
            end
          end

          i=i+1

          verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n", i, tim1,dhim1,him1,ti,dhi,hi)


        end
return (ti,i)

end
