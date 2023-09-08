using YAlpha
using DataFrames
using Plots
using Test


@testset "Potentials" begin
    @test (ħc-197.3269804)<0.1

    h=0.1
    Nmesh=100
    rmesh=0.5*h:h:h*(Nmesh-0.5)
    nu=0.2
    #df_Lambda=read_SkyrmeParam()
	println(df_Lambda)
	ParamIndex=11
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
	qcmesh=10:50:200
	h=0.1
    N_rmesh=100
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
    nu=0.2
	ParamIndex=1

	p1=plot(xlabel="r", ylabel="χ",title="Parameter=$ParamIndex")
	for i=eachindex(qcmesh)
		state=LamAlphaWaveFunc(qcmesh[i],rmesh,nu,ParamIndex)
		p1=plot!(rmesh,state.ψ,label="q=$(state.qc) [MeV/c]")
	end
	savefig(p1,"LamAlphaWaveFunc.pdf")

	p2=plot(xlabel="r", ylabel="χ",title="Spherical bessel function")
	for i=eachindex(qcmesh)
		p2=plot!(rmesh,sin.(qcmesh[i]/ħc*rmesh[:]),label="$(qcmesh[i]) [MeV/c]")
	end
	savefig(p2,"BesselFunction.pdf")

	p3=plot(xlabel="r", ylabel="χ",title="Zero potential")
	zeroarray=zeros(Float64,N_rmesh)
	PSzero=PotSet(ħc^2/(2*940.0)*ones(Float64,N_rmesh),zeroarray,zeroarray,zeroarray)
	for i=eachindex(qcmesh)
		p3=plot!(rmesh,RadWaveFunc(qcmesh[i],940.0,rmesh,PSzero).ψ,label="$(qcmesh[i]) [MeV/c]")
	end
	savefig(p3,"ZeroPotential.pdf")

	qcmesh=50:50:300
	p4=plot(xlabel="q [MeV/c]", ylabel="δ/π")
	δ=zeros(Float64,length(qcmesh))
	for i=eachindex(qcmesh)
		state=LamAlphaWaveFunc(qcmesh[i],rmesh,nu,ParamIndex)
		δ[i]=PhaseShift(state,rmesh)
	end
	p4=plot!(qcmesh,δ)
	savefig(p4,"PhaseShift.pdf")
end

@testset "CoorelationFunction" begin
	qcmesh=10.0:10.0:300.0
	h=0.1
    N_rmesh=1000
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
    nu=0.27
	Pid=2
	R=[1.0,3.0,5.0]
	df_Lambda=read_SkyrmeParam()

	p1=plot(xlabel="q [MeV/c]", ylabel="C(q)", title="ν=$(nu)"
	, left_margin=12Plots.mm, bottom_margin=12Plots.mm)
	for withmom in [true,false]
		for r in R
			C=CoorelationFunction(qcmesh,rmesh,nu,Pid,r,withmom=withmom)
			p1=plot!(qcmesh,C,label="$(df_Lambda[Pid,"ParameterName"]) R= $r fm")
		end
	end
	savefig(p1,"CoorelationFunction.pdf")
end

@testset "LamAlphaBoundState" begin
	h=0.1
    N_rmesh=500
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
    nu=0.27
	ParamIndex=2
	df_Lambda=read_SkyrmeParam()

	E,ψ=Calc_LamAlphaBoundState(rmesh,nu,ParamIndex,withmom=true)

	p1=plot(xlabel="r [fm]", ylabel="ψ", title="B.E. = $(E) MeV"
	, left_margin=12Plots.mm, bottom_margin=12Plots.mm)
	p1=plot!(rmesh,ψ,label="$(df_Lambda[ParamIndex,"ParameterName"])")
	savefig(p1,"BoundState.pdf")

	p2=plot(xlabel="E []", ylabel="ψ", title="B.E. = $(E) MeV"
	, left_margin=12Plots.mm, bottom_margin=12Plots.mm)
	p1=plot!(rmesh,ψ,label="$(df_Lambda[ParamIndex,"ParameterName"])")
	savefig(p1,"BoundState.pdf")
end


@testset "Output" begin
	h=0.1
    N_rmesh=500
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
    nu=0.27
	ParamIndex=[1,2,3,4]

	Output_BoundState(rmesh,nu,ParamIndex)
end