function zoom_Cub(h :: C2LineFunction,
                t₀ :: Float64,
                t₁ :: Float64;
                c₁ :: Float64=0.01,
                tol :: Float64=1e-7,
                ϵ :: Float64=1e-10,
                maxiter :: Int=100,
                verbose :: Bool=false)

        #println("on est dans zoom")

        if obj(h,t₀)<obj(h,t₁)
          tl=t₀
          th=t₁
        else
          tl=t₁
          th=t₀
        end

        #println("t₀=",t₀," t₁=",t₁)

        if t₀<t₁
          tmin=t₀
          tmax=t₁
        else
          tmin=t₁
          tmax=t₀
        end
        #println("tmin=",tmin," tmax=",tmax)

        γ=0.8

        h₀=obj(h,t₀)
        dh₀=grad(h,t₀)

        hl=obj(h,tl)
        dhl=grad(h,tl)
        hh=obj(h,th)
        dhh=grad(h,th)

        t=t₁
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

        # s=t-tqnp
        # y=dhi-dhim1
        #
        # dC=-dhi*s/y
        #
        # tiB=(tim1+t)/2
        # tiC=t+dC
        # i=0
        #
        # #println("tiB=",tiB," tiC=",tiC)
        #
        # if (tiC>tmin) || (tiC<tmax)
        #   #println("on choisit tiC")
        #   ti=tiC
        #   hi=obj(h,tiC)
        #   dhi=grad(h,tiC)
        #   print_with_color(:green,"S")
        # else
        #   #println("on choisit tiB")
        #   ti=tiB
        #   hi=obj(h,tiB)
        #   dhi=grad(h,tiB)
        #   print_with_color(:green,"B")
        # end
        #println("on a les paramètre de zoom")



        verbose && @printf(" iter      tₗ       tₕ        tᵢ        hₗ        hₕ        hᵢ        dhᵢ\n")
        verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)


        while ((i<50) & (abs(dhi)>tol)) || (i==0)
          if i==0
            println("on est dans le while de zoom")
          end

          if (i>0) & ((hi>h₀+c₁*ti*dh₀) || (hl<=hi))
            # println("if de zoomSec")
            # if (hi>h₀+c₁*ti*dh₀)
            #  println("(hi>h₀+c₁*ti*dh₀)")
            # else
            #  println("(hl<=hi)")
            # end
            th=ti
            #tlast=th
          else
            #println("else zoomSec")
            if (abs(dhi)<ϵ) & (i>0)
              # println(abs(dhi),"<",-0.99*dh₀)
              # println("1er if dans le else zoom")
              topt=ti
              iter=i
              return (topt,iter)
            elseif dhi*(th-tl)>=0
              #println("elseif dans le else de zoom")
              th=tl
            end
            tl=ti
            #tlast=ti
          end


          hl=obj(h,tl)
          dhl=grad(h,tl)
          hh=obj(h,th)
          dhh=grad(h,th)

          # him1=hi
          # dhim1=dhi

          # hlast=obj(h,tlast)
          # dhlast=grad(h,tlast)

          s=ti-tqnp
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

          #println("s=",s," y=",y," Γ=,",Γ," yt=",yt)


          tiB=(tim1+ti)/2
          tiC=ti+dC

          him1=hi
          dhim1=dhi

          #println("tiB=",tiB," tiC=",tiC)

          if ((tmin<tiC) & (tiC<tmax)) & (obj(h,tiC)<hl)
            #println("on choisit tiC")
            ti=tiC
            hi=obj(h,tiC)
            dhi=grad(h,tiC)
            verbose && print_with_color(:green,"C")
          else
            #println("on choisit tiB")
            ti=tiB
            hi=obj(h,tiB)
            dhi=grad(h,tiB)
            verbose && print_with_color(:green,"B")
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

          # him1=hi
          # dhim1=dhi


          i=i+1
          #print_with_color(:red,"i=",string(i))
          verbose && @printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n", i,tl,th,ti,hl,hh,hi,dhi)
          #println("th-ti=",th-ti)
          # if th==tl || abs(th-tl)<1e-15
          #   print_with_color(:red,"BREAK!")
          #   break
          # end

        end
        topt=ti
        iter=i


        return (topt,iter)

end
