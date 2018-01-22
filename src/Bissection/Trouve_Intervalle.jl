export trouve_intervalle
function trouve_intervalle(h :: AbstractNLPModel,
                           t₀ :: Float64,
                           inc0 :: Float64;
                           verbose :: Bool = false)

        iter = 1
        h₀ = obj(h, [t₀])[1]
        g₀ = grad(h, [t₀])[1]
        sd = -sign(g₀)
        inc = inc0
        t₁ = t₀ + sd * inc

        h₁ = obj(h, [t₁])[1]
        g₁ = grad(h, [t₁])[1]

        verbose &&
            @printf("iter t₀        g₀        h₀         t₁        ")
        verbose &&    
            @printf("g₀        h₀\n")
        verbose &&
            @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n",
                iter,  t₀, g₀, h₀, t₁, g₀, h₀)

        while (g₁ * sd < 0.0) && (h₁ < h₀) && (iter < 30)
          inc = inc * 4
          t₀ = t₁
          h₀ = h₁
          g₀ = g₁
          t₁ = t₀ + sd * inc
          h₁ = obj(h, [t₁])[1]
          g₁ = grad(h, [t₁])[1]
          verbose && @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n",
                              iter,  t₀, g₀, h₀, t₁, g₀, h₀)
          iter = iter + 1
        end

        while (g₁ * sd < 0.0) & (iter < 30)
          tₘ = (t₁ + t₀) / 2
          hₘ = obj(h, [tₘ])[1]
          gₘ = grad(h, [tₘ])[1]
          if gₘ * sd>0
            t₁ = tₘ
            g₁ = gₘ
            h₁ = hₘ
          else
            if hₘ < h₀
              t₀ = tₘ
              h₀ = hₘ
              g₀ = gₘ
            else
              t₁ = tₘ
              h₁ = hₘ
              g₁ = gₘ
            end
          end
          verbose &&
            @printf("%4d %7.2e %7.2e  %7.2e  %7.2e  %7.2e  %7.2e \n",
                    iter,  t₀, g₀, h₀, t₁, g₀, h₀)
          iter = iter + 1
        end

        a = min(t₀, t₁)
        b = max(t₀, t₁)

        return (a, b)

end
