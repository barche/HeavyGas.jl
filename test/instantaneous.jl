using HeavyGas
using Interpolations

const dt = 0.01

m = GasData(6.4, 10, 0.003, 20+273.15, 44.10, HeavyGas.D)
@show t_box, x_cen, r, h, vfrac_box, u_box = instantaneous(7., 4.12, 1., dt, m)

# Compute observer positions
tx = x_cen-r
lx = x_cen+r
obs_dx = r/11
obs_xs = tx+obs_dx:obs_dx:lx-obs_dx
vfracs = Vector{Interpolations.GriddedInterpolation{Float64,1,Float64,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1}},0}}(length(obs_xs))

for (i,x) in enumerate(obs_xs)
  Beff0 = sqrt(r^2 - (x_cen-x)^2)
  x, Meff, Beff, Sy, ypol, b, Heff, ueff = dispersion(0.1, 640., Beff0, m, H0 = h, ypol0 = vfrac_box, u0 = u_box)
  vfracs[i] = interpolate((x,), ypol, Gridded(Linear()))
end
#x, Meff, Beff, Sy, ypol, b, Heff, ueff = dispersion(0.1, 580., 46.146, m, H0 = 1.0419, ypol0 = 9.0920e-2, u0 = 2.64)



for t in 0.:20.:100.
  println("\n t = $(t+28)\n")
  for (i,xs) in enumerate(obs_xs)
    x = observer_position(xs, t, r, h, m)
    println(round(x), " ", round(100*vfracs[i][x-xs],2))
    #println(x-xs, " ", vfracs[i][x-xs])
  end
end
