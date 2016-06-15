using Base.Test
using HeavyGas

# Read a HGSYSTEM output file
function readhgfile(filename)
    x = Float64[]
    y = Float64[]
    f = open(filename)
    for ln in eachline(f)
        sl = split(ln)
        if !isempty(sl) && isa(parse(sl[1]), Number)
            push!(x,Float64(parse(sl[1])))
            push!(y,Float64(parse(sl[2])))
        end
    end
    close(f)
    return x, y
end

m = GasData(3, 10, 0.03, 15+273.15, 40.53, HeavyGas.D)
x, Meff, Beff, Sy, ypol, b, Heff, ueff = dispersion(0.1, 497.4, 2.5, m, E0 = 1.)

hgpath = joinpath(dirname(@__FILE__), "hgsystem") # Base path to HGSYSTEM reference files
x_hg, vfrac_hg = readhgfile(joinpath(hgpath, "HTAG01.PS1"))
x_hg, Sy_hg = readhgfile(joinpath(hgpath, "HTAG01.PS2"))
x_hg, Beff_hg = readhgfile(joinpath(hgpath, "HTAG01.PS3"))
x_hg, Sz_hg = readhgfile(joinpath(hgpath, "HTAG01.PS4"))
x_hg, b_hg = readhgfile(joinpath(hgpath, "HTAG01.PS6"))
x_hg, Heff_hg = readhgfile(joinpath(hgpath, "HTAG01.PS8"))

x_hg -= 2.5 # Offset of 2.5 m in HG ref data
const Vm_hg = 8314.3*(15+273.15)/101325
const α = 0.22664
const β = 1 + α
Meff_hg = 2*m.uref*gamma((1+α)/β)/(β*m.href^α*Vm_hg)*Beff_hg.*Sz_hg.^(1+α)

@test x_hg ≈ collect(x)
@test maximum(abs(vfrac_hg-(100*ypol))) < 2.
@test maximum(abs(Sy_hg-Sy)) < 0.2
@test maximum(abs(Beff_hg-Beff)) < 0.3
@test maximum(abs(b_hg-b)) < 0.3
@test maximum(abs(Heff_hg-Heff)) < 0.7
@test maximum(abs(Meff_hg-Meff)) < 12.

# using Plots
# pyplot()
# plot(x, Beff, label="Julia")
# plot!(x_hg, Beff_hg, label="HG")
# show()
