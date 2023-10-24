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


function CorrelationFunction(qcmesh,rmesh,nu,ParamIndex::Int,R,df_Lambda; withmom=true, Gauss=false)
	C=zeros(Float64,length(qcmesh))

	for i=eachindex(qcmesh)
		C[i]=CalcCF(qcmesh[i],rmesh,nu,ParamIndex,R,df_Lambda,withmom=withmom,Gauss=Gauss)
	end

	return C
end


function ScatAmp_LowE(qcmesh,a0,reff)
	return (@. (-1/a0 + reff*qcmesh[:]^2/(2*ħc^2) - im*qcmesh[:]/ħc)^(-1))
end


function Func1_Int(x::AbstractFloat)
	f(t)=exp(t^2-x^2)/x
	tmesh=range(0,x,1000)
	return MyLib.IntTrap(tmesh,f.(tmesh))
end

function Func1(x::AbstractArray)
	return Func1_Int.(x)
end


function Func2(x)
	return (@. (1-exp(-x^2))/x)
end

function Func3(x)
	return 1-x/(2*π^0.5)
end

function LLformula_S1(qcmesh,a0,reff,R)
	fq = ScatAmp_LowE(qcmesh,a0,reff)
	x=2*R*qcmesh/ħc
	F1 = Func1(x)
	F2 = Func2(x)
	F3 = Func3(reff/R)

	CF = ones(Float64,length(qcmesh))
	@. CF += abs2(fq)*F3/(2*R^2)
	@. CF += 2*real(fq)*F1/(π^0.5*R)
	@. CF -= imag(fq)*F2/R

	return CF
end

export CorrelationFunction, LLformula_S1

end