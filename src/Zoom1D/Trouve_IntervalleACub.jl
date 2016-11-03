function trouve_intervalleACub(h :: C2LineFunction,
                t₀ :: Float64,
                tₘ :: Float64;
                ϵ :: Float64=1e-10,
                verbose :: Bool=true)

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

        hi=obj(h,t)
        dhi=grad(h,t)
        him1=obj(h,tqnp)
        dhim1=grad(h,tqnp)
        #
        # s=t-tqnp
        # y=dhi-dhim1
        #
        # println("s=",s," y=",y)
        #
        # dS=-dhi*s/y
        #
        # tiB=(tim1+t)/2
        # tiC=t+dS
        #
        # println("tiB=",tiB," tiC=",tiC)
        # i=0
        #
        # if (tiC<t₀) || (tiC>tₘ) || (tiC==NaN)
        #   println("on choisit tiB")
        #   ti=tiB
        #   hi=obj(h,tiB)
        #   dhi=grad(h,tiB)
        #   print_with_color(:green,"B")
        # else
        #   println("on choisit tiC")
        #   ti=tiC
        #   hi=obj(h,tiC)
        #   dhi=grad(h,tiC)
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

        println("on a les paramètres de trouve intervalle")

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
            (topt,iter)=zoom_Cub(h,tim1,ti)
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
            (topt,iter)=zoom_Cub(h,ti,tim1)
            # if obj(h,topt)<obj(h,t₀)
            #   println("h(t*)<h(t₀)")
            # end
            return (topt,iter)
          end

          s=t-tqnp
          y=dhi-dhim1


          α=-s
          z=dhi+dhim1+3*(hi-him1)/α
          discr=z^2-dhi*dhim1
          denom=dhi+dhim1+2*z

          if (discr>0) & (abs(denom)>eps(Float64))
            #si on peut on utilise l'interpolation cubique
            #println("si on peut on utilise l'interpolation cubique")
            w=sqrt(discr)
            dC=-s*(dhi+z+sign(α)*w)/(denom)
          else #on se rabat sur une étape de
            dC=-dhi*s/y
          end


          tiB=(ti+t)/2
          tiC=t+dC

          him1=hi
          dhim1=dhi

          #println("tiB=",tiB," tiC=",tiC)

          if (tiC<t₀) || (tiC>tₘ)
            #println("on choisit tiB")
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            print_with_color(:green,"B")
          else
            #println("on choisit tiC")
            ti=tiC
            hi=obj(h,tiC)
            dhi=grad(h,tiC)
            print_with_color(:green,"C")
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
