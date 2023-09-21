module CorrelationFunc

include("./Scattering.jl")
using .Scattering
include("./LamAlphaPot.jl")
import .LamAlphaPot
include("./constants.jl")
include("./MyLib.jl")
import .MyLib

mutable struct CF
	qc::AbstractArray
	C::Vector{AbstractArray}
end


function LamAlphaWaveFunc(qc,rmesh,nu,ParamIndex,df_Lambda; withmom=true, Gauss=false)
	LamAPot=LamAlphaPot.CalcPotentials(rmesh,nu,ParamIndex,df_Lambda,Gauss=Gauss)
	μ=mΛMeV*mαMeV/(mΛMeV + mαMeV)
	if withmom==false
		LamAPot.h2_2μeff=ħc^2/(2*μ)*ones(Float64,length(rmesh))
		LamAPot.dh2_2μeff=zeros(Float64,length(rmesh))
		LamAPot.ddh2_2μeff=zeros(Float64,length(rmesh))
	end
	PS=Scattering.PotSet(LamAPot.h2_2μeff, LamAPot.dh2_2μeff, LamAPot.ddh2_2μeff, LamAPot.U_local)
	state=Scattering.RadWaveFunc(qc,μ,rmesh,PS)
	return state
end

export LamAlphaWaveFunc

function CalcCF(qc,rmesh,nu,ParamIndex,R::AbstractFloat,df_Lambda; withmom=true, Gauss=false)
	st=LamAlphaWaveFunc(qc,rmesh,nu,ParamIndex,df_Lambda,withmom=withmom,Gauss=Gauss)
	y=@. exp(-rmesh[:]^2/(4*R^2))
	@. y*=(st.ψ[:]^2 - (sin(qc/ħc*rmesh[:]))^2)
	I_source=MyLib.IntTrap(rmesh,y)
	I_source+=rmesh[1]*y[1]*0.5
	return 1+I_source/(2*π^0.5*R^3*qc^2/ħc^2)
end

function CoorelationFunction(qcmesh,rmesh,nu,ParamIndex::Int,R,df_Lambda; withmom=true, Gauss=false)
	C=zeros(Float64,length(qcmesh))

	for i=eachindex(qcmesh)
		C[i]=CalcCF(qcmesh[i],rmesh,nu,ParamIndex,R,df_Lambda,withmom=withmom,Gauss=Gauss)
	end

	return C
end

export CoorelationFunction

end