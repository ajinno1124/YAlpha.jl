module CoorelationFunc

include("./Scattering.jl")
import .Scattering
include("./LamAlphaPot.jl")
import .LamAlphaPot
include("./constants.jl")

mutable struct CF
	q::AbstractArray
	C::Vector{AbstractArray}
end

function LamAlphaWaveFunc(q,rmesh,nu,ParamIndex)
	LamAPot=LamAlphaPot.CalcPotentials(rmesh,nu,ParamIndex)
	PS=Scattering.PotSet(LamAPot.h2_2μeff, LamAPot.dh2_2μeff, LamAPot.ddh2_2μeff, LamAPot.U_local)
	μ=mΛMeV*mαMeV/(mΛMeV + mαMeV)
	state=Scattering.RadWaveFunc(q,μ,rmesh,PS)
	return state
end

export LamAlphaWaveFunc

function CalcCF(q,rmesh,nu,ParamIndex)
	χ=LamAlphaWaveFunc(q,rmesh,nu,ParamIndex)

	return C
end

function CoorelationFunction(nu,ParamIndex::Int)
	qmesh=0.1:0.1:5
	h=0.1
	N_rmesh=200
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
	C=zeros(Float64,length(qmesh))

	for i=eachindex(qmesh)
		C[i]=CalcCF(qmesh[i],rmesh,nu,ParamIndex)
	end

	return CF(qmesh,C)

end

end