using YAlpha
using DataFrames
using Plots
using Test

@testset "Potentials" begin
    @test (ħc-197.3269804)<0.1

    h=0.1
    Nmesh=1000
    rmesh=0.5*h:h:h*(Nmesh-0.5)
    nu=0.2
    df_Lambda=read_SkyrmeParam("Lambda Parameters.dat")
    aL=getaL(df_Lambda,1)
    @test aL[1]==-500.89

    DensePots=CalcPotentials(rmesh,nu,aL)
    plot(xlabel="r (fm)", title="LYIV",xlim=(0,3))
    plot!(rmesh,DensePots.ρ,label="ρ")
    plot!(rmesh,DensePots.τ,label="τ")
    plot!(rmesh,DensePots.Lapρ,label="Lapρ")
    plot!(rmesh,DensePots.Ulocal,label="U_local")
    savefig("Densities.pdf")

    plot(xlabel="r (fm)", title="LYIV",xlim=(0,3))
    plot!(rmesh,DensePots.h2_2μeff,label="h2_2μeff")
    plot!(rmesh,DensePots.dh2_2μeff,label="dh2_2μeff")
    plot!(rmesh,DensePots.ddh2_2μeff,label="ddh2_2μeff")
    savefig("EffectiveMass.pdf")
end

