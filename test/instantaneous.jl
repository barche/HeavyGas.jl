using HeavyGas
using Interpolations

const dt = 0.01

m = GasData(6.4, 10, 0.003, 20+273.15, 44.10, HeavyGas.D)

x, vfrac = instantaneous(7., 4.12, 1., dt, 28.:20.:128., m, 2)
println(x,"\n",vfrac)
