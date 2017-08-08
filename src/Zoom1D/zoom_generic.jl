export zoom_generic

function zoom_generic(h :: LineModel,
                  	  t₀ :: Float64,
	                  t₁ :: Float64;
					  direction :: String = "Nwt",
	                  tol :: Float64 = 1e-7,
	                  ϵ :: Float64 = 1e-10,
	                  maxiter :: Int = 100,
	                  verbose :: Bool = false,
					  kwargs...)
    if obj(h, t₀) < obj(h, t₁)
    tl = t₀; th = t₁
    else
    tl = t₁; th = t₀
    end

    hl = obj(h, tl); dhl = grad(h, tl)
    hh = obj(h, th); dhh = grad(h, th)
    h₀ = obj(h, t₀); dh₀ = grad(h, t₀)

    γ = 0.8; t = t₁; tp = t₀; tqnp = t₀; hk = 0; hkm1 = 0; gkm1 = 0; hp = 0
    i = 0

    hk = obj(h, t); gk = grad(h, t)
    hkm1 = obj(h, tqnp); gkm1 = grad(h, tqnp)

    verbose &&
        @printf(" iter        tₗ        tₕ         t        hₗ        hₕ")
        @printf("         hk         gk\n")
    verbose &&
        @printf(" %7.2e %7.2e  %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e\n",
                i, tl, th, t, hl, hh, hk, gk)


    while ((i < maxiter) && (abs(gk) > tol)) || i == 0

    if (hk > h₀) || (hl <= hk)
      th = t; tlast = th; hlast = hh
      dhlast = dhh; hh = hk; dhh = gk
    else
      if (abs(gk) < ϵ)
        topt = t; iter = i
        return (topt, iter)
      elseif gk * (th - tl) >= 0
        th = tl; hh = hl; dhh = dhl
      end
      tl = t; tlast = tl; hlast = hl
      dhlast = dhl; hl = hk; dhl = gk
    end
    if direction == "Nwt"
        kₖ = hess(h, t)
        dN = -gk/kₖ #direction de Newton
    elseif direction == "Sec"
        s = t - tqnp; y = gk - gkm1
        dN =-gk * s / y
    elseif direction == "SecA"
        s = t - tqnp; y = gk - gkm1
        Γ = 3 * (gk + gkm1) * s - 6 * (hk - hkm1)
        if y * s + Γ < eps(Float64) * (s^2)
          yt = y
        else
          yt = y + Γ / s
        end
        dN = -gk * s / yt
    elseif direction == "Cub"
        s = t - tqnp; y = gk - gkm1; α = -s
        z = gk + gkm1 + 3 * (hk - hkm1) / α
        discr = z^2 - gk * gkm1; denom = gk + gkm1 + 2 * z
        if (discr > 0.0) && (abs(denom) > eps(Float64))
          w = sqrt(discr)
          dN = -s * (gk + z + sign(α) * w) / (denom)
        else
          dN=-gk*s/y
        end
    end

    if ((tp - t) * dN > 0) & (dN/(tp - t) < γ)
      tplus = t + dN
      hplus = obj(h, tplus); gplus = grad(h, tplus)
      verbose && println("N")
    else
      tplus = (t + tp) / 2
      hplus = obj(h, tplus); gplus = grad(h, tplus)
      verbose && println("B")
    end

    if t > tp
      if gplus < 0.0
        tp = t; tqnp = t; t = tplus
      else
        tqnp = t; t = tplus
      end
    else
      if gplus > 0.0
        tp = t; tqnp = t; t = tplus
      else
        tqnp = t; t = tplus
      end
    end

    #mise à jour des valeurs
    hkm1 = hk; gkm1 = gk; hk = hplus; gk = gplus
    i = i + 1
    #println("on continue d'itérer")
    verbose &&
        @printf(" %7.2e %7.2e  %7.2e  %7.2e %7.2e %7.2e %7.2e %7.2e\n",
                i, tl, th, t, hl, hh, hk, gk)

    end
    topt = t; iter = i

    return (topt, iter)
end #function
