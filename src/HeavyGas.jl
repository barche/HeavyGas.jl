module HeavyGas

using Interpolations
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

  function GasData(uref, href, z0, Tg, mdp, class::StabilityClass)
    @assert class == D
    # TODO: Only for class D for now
    δ600 = 0.08
    λ = 1e5
    Ψ(z) = -6.9*z/λ
    uτ = κ*uref/(log((href+z0)/z0) - Ψ(href))
    zvals = collect(0.:0.25:40.)
    fit = curve_fit((z,p) -> uref*(z/href).^p[1], zvals, Float64[uτ/κ*(log((z+z0)/z0) - Ψ(z)) for z in zvals], Float64[1/(1+10*z/href) for z in zvals], [1.])

    const P = 101325. # pressure (Pa)
    const R = 8314.3 # gas constant (J/(kmol*K))
    const ma = 28.95 # Molecular mass of air (kg / kmol)

    ρa = P*ma/(R*Tg)

    new(uref, href, z0, mdp, δ600, uτ, fit.param[end], 1+fit.param[end], ma, ρa, ma/ρa, λ)
  end
end

ρ(d::GasData, z::Float64) = d.ρa*(1-d.uτ^2 / (κ^2*g*d.λ) * (log((z+d.z0)/d.z0) + 9.2*z/d.λ))

function dispersion(dx::Float64, xend::Float64, Beff0::Float64, d::GasData; E0 = 0., ypol0 = 0., H0 = 0., u0 = 0.)
  const ut = d.uτ
  const CE = 1.15
  const tav = 600 # Averaging time
  const δ = d.δ600*(tav/600)
  const γ = 0.0001

  const ca0 = d.mdp*ypol0/d.Vm
  const Sz0 = H0 / (1/d.β*gamma(1/d.β))
  const E = E0 != 0 ? E0 : (2*d.uref*gamma((1+d.α)/d.β) / (d.β*d.href^d.α) * ca0 * Beff0) * Sz0^(1+d.α)
  const Meff0 = (E / d.mdp) / (ypol0 == 0 ? 1 : ypol0)

  xrange = dx:dx:xend

  Meff_arr = Vector{Float64}(length(xrange)+1)
  Beff_arr = Vector{Float64}(length(xrange)+1)
  Sy_arr = Vector{Float64}(length(xrange)+1)
  b_arr = Vector{Float64}(length(xrange)+1)
  Heff_arr = Vector{Float64}(length(xrange)+1)
  ueff_arr = Vector{Float64}(length(xrange)+1)

  Meff_arr[1] = Meff0
  Beff_arr[1] = Beff0
  Sy_arr[1] = 0
  b_arr[1] = Beff0

  gaussian = false
  gravity_collapsed = false
  xv = 0.

  for (i,x) in enumerate(xrange)
    Meff = Meff_arr[i]
    Beff = Beff_arr[i]
    Sy = Sy_arr[i]
    ypol = E/(d.mdp*Meff) # polutant molar fraction
    ca = d.mdp*ypol/d.Vm
    ρm = ca + (1-ypol)*d.ma / d.Vm
    Sz = (E / (2*d.uref*gamma((1+d.α)/d.β) / (d.β*d.href^d.α) * ca * Beff))^(1/(1+d.α))
    Heff = 1/d.β*gamma(1/d.β)*Sz
    Heff_arr[i] = Heff
    ueff = gamma((1+d.α)/d.β)/gamma(1/d.β)*d.uref*(Sz/d.href)^d.α
    ueff_arr[i] = ueff
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

    if !gravity_collapsed && ((Beff/Heff) / sqrt(Ri*(1+0.8*(Ris))) >= 8/(3*κ))
      gravity_collapsed = true
    end

    b_arr[i+1] = b
    if gravity_collapsed
      Meff_arr[i+1] = Meff/(2*Beff) + dx*(ue/V0)
    else
      Meff_arr[i+1] = Meff + dx*(2*Beff*ue/V0)
    end

    ke(S::Float64) = 2*δ^2/γ*(1/S-1/S^2)

    Beff_pr = 0.
    if b > 0 && Ris > 0
      if gravity_collapsed
        S = 1 + sqrt(1+8*(δ/(γ*Sy))^2)
        Beff_arr[i+1] = sqrt(π*ke(S)*dx + Beff^2)
        Meff_arr[i+1] *= 2*Beff_arr[i+1]
      else
        Beff_pr = CE * d.uτ*sqrt(Ri)/ueff
        Beff_arr[i+1] = Beff + dx*Beff_pr
      end
    else
      Beff_arr[i+1] = Beff
    end

    S = 1 + sqrt(1+2*π*(δ/(γ*Beff))^2)
    Sy_arr[i+1] = b > 0 ? sqrt(4*ke(S)*dx + Sy^2) : sqrt(2)*δ*(x+xv)/sqrt(1+γ*(x+xv))
  end

  Heff_arr[end] = Heff_arr[end-1]
  ueff_arr[end] = ueff_arr[end-1]

  return 0.:dx:xend, Meff_arr, Beff_arr, Sy_arr, E./(d.mdp*Meff_arr), b_arr, Heff_arr, ueff_arr
end

# Compute initial gravitational spread using a box model (Hegabox)
function boxmodel(r0::Number, h0::Number, c0::Number, dt::Float64, d::GasData)
  r::Float64 = r0
  h::Float64 = h0
  V = π*r^2*h
  t = 0.
  ub = 0.
  x = 0.

  const Mdp = V*c0/d.Vm
  const αe = 0.85
  const ts = 1.75#3.6#0.83

  function ρ(V::Float64)
    c = Mdp*d.Vm / V
    ca = d.mdp*c/d.Vm
    return ca + (1-c)*d.ma / d.Vm
  end

  Ris_i = 1e30

  const delay_t = ts / sqrt((g*(ρ(V)-d.ρa)/d.ρa) / V^(1/3))

  Ris_i = 100.
  while Ris_i > 10
    m = ρ(V)*V
    mv = ub*m
    gp = g*(ρ(V) - d.ρa)/d.ρa
    uf = 1.15*sqrt(gp*h)
    r += dt * uf

    ue = αe*uf
    Ris = gp*h/d.uτ^2
    int_u = d.uτ/κ*(h*log((h + d.z0)/d.z0) - h - d.z0*log(d.z0) + d.z0*log(h + d.z0))
    int_u / h
    ua = (0.8 + 0.2/(1+Ris))/h * int_u
    usi = κ*ub/(log(h/d.z0)-1)
    Ris_i = usi == 0 ? 1e30 : gp*h/usi^2
    ϕ = sqrt(1+0.8*Ris_i)
    ut = κ*usi/ϕ
    dV = t > delay_t ? dt*(2*π*r*h*ue + π*r^2*ut) : 0.
    V += dV
    ub = (mv + dV*d.ρa*ua) / (m + dV*d.ρa)
    h = V / (π*r^2)
    t += dt
    x += ub*dt
  end
  return t, x, r, h, Mdp*d.Vm / V, ub
end

function distribute_observers(x::Float64, r::Float64, n::Int)
  const dx = r/1e6
  return linspace(x-r+dx, x+r-dx, 2*n+1)
end

"Data linked to observers"
type ObserverDataBase{InterpolationType}
  xstart::AbstractArray{Float64,1}
  vfrac::Vector{InterpolationType}
  b::Vector{InterpolationType}
  Beff::Vector{InterpolationType}

  function ObserverDataBase(x::Float64, r::Float64, n::Int)
    const xrange = distribute_observers(x,r,n)
    const nx = length(xrange)
    new(xrange, Vector{InterpolationType}(nx), Vector{InterpolationType}(nx), Vector{InterpolationType}(nx))
  end
end

typealias ObserverData ObserverDataBase{Interpolations.GriddedInterpolation{Float64,1,Float64,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1}},0}}

function observer_position(xs::Float64, t::Float64, Rm::Float64, Heff::Float64, d::GasData)
  const Sz = d.β/gamma(1/d.β)*Heff
  const C0 = d.uref*gamma((1+d.α)/d.β)/gamma(1/d.β) * (Sz/d.href)^d.α * ((0.5*π*Rm + Rm)/d.href)^(-d.α/(1+d.α))
  return xs + d.href*(C0/(d.href*(1+d.α))*(t))^(1+d.α)
end

function correct_cloud_shape(observer_positions::Vector{Float64})
  Lc = observer_positions[end] - observer_positions[1]
end

function instantaneous(r0::Number, h0::Number, c0::Number, dt::Float64, output_times::AbstractArray{Float64,1}, d::GasData, n_obs = 10, dx = 0.1)
  t_box, x_box, r_box, h_box, vfrac_box, u_box = boxmodel(r0, h0, c0, dt, d)
  od = ObserverData(x_box, r_box, n_obs)

  const tmax = maximum(output_times)
  xmax = 1.1*observer_position(od.xstart[end], tmax, r_box, h_box, d)

  Beff0 = zeros(od.xstart)

  for (i,x) in enumerate(od.xstart)
    Beff0[i] = sqrt(r_box^2 - (x_box-x)^2)
    x, Meff, Beff, Sy, ypol, b, Heff, ueff = dispersion(dx, xmax, Beff0[i], d, H0 = h_box, ypol0 = vfrac_box, u0 = u_box)
    od.vfrac[i] = interpolate((x,), ypol, Gridded(Linear()))
    od.b[i] = interpolate((x,), b, Gridded(Linear()))
    od.Beff[i] = interpolate((x,), Beff, Gridded(Linear()))
  end

  const Lc = od.xstart[end] - od.xstart[1]
  const Lcs = Lc / r_box
  Wc = zeros(od.xstart)
  Ss = zeros(od.xstart)
  Beff_c = zeros(od.xstart)
  x_corr = zeros(od.xstart)

  for t in output_times
    x_orig = Float64[observer_position(xs, t-t_box, r_box, h_box, d) for xs in od.xstart]
    println("\n t = $(t)\n")
    for (i,xs) in enumerate(od.xstart)
      x = x_orig[i]
      const xo = x - xs
      Wc[i] = od.b[i][xo] > 0 ? od.Beff[i][xo] / Beff0[i] : 1.
      Ss[i] = -0.5*(1+Lcs) + 0.5*sqrt((1+Lcs)^2 + 4*Lcs*(Wc[i]-1))
      Beff_c[i] = od.Beff[i][xo] * (1+Ss[i]) / Wc[i]
      println(round(x), " ", round(100*od.vfrac[i][xo],2), " ", od.b[i][xo], " ", od.Beff[i][xo])
    end

    const xend = endof(od.xstart)
    const xmid = endof(od.xstart)÷2+1
    x_corr[xmid] = x_orig[xmid]
    for i in xmid-1:-1:1
      x_corr[i] = x_corr[i+1] + (x_orig[i] - x_orig[i+1])*(Lcs + Ss[i])/Lcs
    end
    for i in xmid+1:xend
      x_corr[i] = x_corr[i-1] - (x_orig[i-1] - x_orig[i])*(Lcs + Ss[i])/Lcs
    end

    println("corrected cloud:")

    for (i,xs) in enumerate(od.xstart)
      x = x_corr[i]
      const xo = x - xs
      println(round(x), " ", round(100*od.vfrac[i][xo],2), " ", od.b[i][xo], " ", od.Beff[i][xo], " ", Beff_c[i])
    end

  end
end

export GasData, dispersion, boxmodel, ObserverData, instantaneous

end # module
