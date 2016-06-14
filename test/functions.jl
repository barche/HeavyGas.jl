using Base.Test
using HeavyGas
using Interpolations

for n in [1,2,3,4,10,11,100]
  for x in [0., 1., π]
    for r in [1., 11., sqrt(2)]
      @test length(HeavyGas.distribute_observers(x, r, n)) == 2*n+1
      @test HeavyGas.distribute_observers(x, r, n)[end÷2+1] ≈ x
    end
  end
end

od = ObserverData(0., 10., 10)
@test length(od.xstart) == 21
@test od.xstart[end÷2+1] ≈ 0.
