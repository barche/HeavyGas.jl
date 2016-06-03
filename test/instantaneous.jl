using HeavyGas

m = GasData(6.4, 10, 0.003, 20+273.15, 44.10, 32, HeavyGas.D)
instantaneous(7., 4.12, 1., 0.001, m)
