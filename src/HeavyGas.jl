module HeavyGas

using LsqFit

# Pasquill stability classes
@enum StabilityClass A B C D E F

const V0 = 22.4 # m³/kmol
const g = 9.81
const κ = 0.41

immutable GasData
  "Reference wind speed (m/s)"
  uref::Float64

  "Altitude for reference wind speed (m)"
  href::Float64

  "Surface rougness (m)"
  z0::Float64

  "Molecular mass of dry pollutant (kg/kmol)"
  mdp::Float64

  "Source strength (kg of dp/s)"
  E::Float64

  δ600::Float64

  "Friction velocity (m/s)"
  uτ::Float64

  "Power law fit parameter"
  α::Float64

  "Concentration profile parameter"
  β::Float64

  "Molecular mass of air"
  ma::Float64

  "Air density at ground level"
  ρa::Float64

  "Molecular volume"
  Vm::Float64

  "Monin-Obukhov length (m)"
  λ::Float64

  function GasData(uref, href, z0, Tg, mdp, E, class::StabilityClass)
    @assert class == D
    # TODO: Only for class D for now
    δ600 = 0.08
    λ = 1e5
    Ψ(z) = -6.9*z/λ
    @show uτ = κ*uref/(log((href+z0)/z0) - Ψ(href))
    zvals = collect(0.:0.25:40.)
    fit = curve_fit((z,p) -> uref*(z/href).^p[1], zvals, Float64[uτ/κ*(log((z+z0)/z0) - Ψ(z)) for z in zvals], Float64[1/(1+10*z/href) for z in zvals], [1.])
    #fit.param[end] = 0.22664

    const P = 101325. # pressure (Pa)
    const R = 8314.3 # gas constant (J/(kmol*K))
    const ma = 28.95 # Molecular mass of air (kg / kmol)

    new(uref, href, z0, mdp, E, δ600, uτ, fit.param[end], 1+fit.param[end], ma, P*ma/(R*Tg), R*Tg/P, λ)
  end
end

ρ(d::GasData, z::Float64) = d.ρa*(1-d.uτ^2 / (κ^2*g*d.λ) * (log((z+d.z0)/d.z0) + 9.2*z/d.λ))

function dispersion(dx::Float64, xend::Float64, Meff0::Float64, Beff0::Float64, d::GasData)
  const ut = d.uτ
  const CE = 1.15
  const tav = 600 # Averaging time
  const δ = d.δ600*(tav/600)
  const γ = 0.0001

  xrange = dx:dx:xend

  Meff_arr = Vector{Float64}(length(xrange)+1)
  Beff_arr = Vector{Float64}(length(xrange)+1)
  Sy_arr = Vector{Float64}(length(xrange)+1)
  b_arr = Vector{Float64}(length(xrange)+1)
  Heff_arr = Vector{Float64}(length(xrange)+1)

  Meff_arr[1] = Meff0
  Beff_arr[1] = Beff0
  Sy_arr[1] = 0
  b_arr[1] = Beff0

  println("Running with α = $(d.α)")

  gaussian = false
  gravity_collapsed = false
  xv = 0.

  for (i,x) in enumerate(xrange)
    Meff = Meff_arr[i]
    Beff = Beff_arr[i]
    Sy = Sy_arr[i]
    ypol = d.E/(d.mdp*Meff) # polutant molar fraction
    ca = d.mdp*ypol/d.Vm
    ρm = ca + (1-ypol)*d.ma / d.Vm
    Sz = (d.E / (2*d.uref*gamma((1+d.α)/d.β) / (d.β*d.href^d.α) * ca * Beff))^(1/(1+d.α))
    Heff = 1/d.β*gamma(1/d.β)*Sz
    Heff_arr[i] = Heff
    ueff = gamma((1+d.α)/d.β)/gamma(1/d.β)*d.uref*(Sz/d.href)^d.α
    ρh = ρ(d, Heff)
    Ris = g*Heff/ut^2*((ρm - ρh)/d.ρa)
    Ri = g*Heff/d.uτ^2*(1-d.ρa/ρm)
    ϕ = Ris > 0 ? sqrt(1+0.8*Ris)/(1+d.α) : 1/(sqrt(1-0.6*Ris)*(1+d.α))
    ue = κ*ut/ϕ
    b = Beff - 0.5*sqrt(π)*Sy
    if b <= 0 && !gaussian
      println("Profile becomes gaussian at x = $x")
      gaussian = true
      xv = (Sy^2*γ - 4*δ^2*x + sqrt(Sy^2*(Sy^2*γ^2 + 8*δ^2)))/(4*δ^2)
    end
    if gaussian
      b = 0.
    end

    b_arr[i+1] = b
    Meff_arr[i+1] = Meff + dx*(2*Beff*ue/V0)

    if !gravity_collapsed && ((Beff/Heff) / sqrt(Ri*(1+0.8*Ris)) >= 8/(3*κ))
      gravity_collapsed = true
    end
    Beff_pr = 0.
    if b > 0 && Ris > 0
      Beff_pr = gravity_collapsed ? d.uτ/ueff * Ri*sqrt(1+0.8*Ris)/(3*κ*5) * (Heff/Beff) : CE * d.uτ*sqrt(Ri)/ueff
    end
    Beff_arr[i+1] = Beff + dx*Beff_pr

    σy = sqrt(2/pi)*Beff
    xe = σy*(γ*σy + sqrt(4*δ^2 + γ^2*σy^2))/(2*δ^2)
    ke = σy*(δ/sqrt(1+γ*xe) - 0.5*δ*xe*γ*(1+γ*xe)^(-3/2))
    Sy_arr[i+1] = b > 0 ? sqrt(4*ke*dx + Sy^2) : sqrt(2)*δ*(x+xv)/sqrt(1+γ*(x+xv))
  end

  Heff_arr[end] = Heff_arr[end-1]

  return 0.:dx:xend, Meff_arr, Beff_arr, Sy_arr, d.E./(d.mdp*Meff_arr), b_arr, Heff_arr
end

export GasData, dispersion

end # module
