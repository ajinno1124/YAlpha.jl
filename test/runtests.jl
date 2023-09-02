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
    df_Lambda=read_SkyrmeParam()
	ParamIndex=1
    aL=getaL(ParamIndex)
    @test aL[1]==-500.89

    DensePots=CalcPotentials(rmesh,nu,ParamIndex)
    plot(xlabel="r (fm)", title="LYIV",xlim=(0,3))
    plot!(rmesh,DensePots.ρ,label="ρ")
    plot!(rmesh,DensePots.τ,label="τ")
    plot!(rmesh,DensePots.Lapρ,label="Lapρ")
    plot!(rmesh,DensePots.U_local,label="U_local")
    savefig("Densities.pdf")

    plot(xlabel="r (fm)", title="LYIV",xlim=(0,3))
    plot!(rmesh,DensePots.h2_2μeff,label="h2_2μeff")
    plot!(rmesh,DensePots.dh2_2μeff,label="dh2_2μeff")
    plot!(rmesh,DensePots.ddh2_2μeff,label="ddh2_2μeff")
    savefig("EffectiveMass.pdf")
end

@testset "Scattering" begin
	qcmesh=50:50:200
	h=0.1
    N_rmesh=300
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
    nu=0.2
	ParamIndex=1
	plot(xlabel="r", ylabel="χ",title="Parameter=$ParamIndex")
	for i=eachindex(qcmesh)
		state=LamAlphaWaveFunc(qcmesh[i],rmesh,nu,ParamIndex)
		plot!(rmesh,state.ψ,label="q=$(state.qc) [MeV/c]")
	end
	savefig("YAlphaPotential.pdf")
end

