module CoorelationFunc

include("./Scattering.jl")
import .Scattering
include("./LambdaAlphaPot.jl")
import .LambdaAlphaPot

mutable struct CF
	q::AbstructArray
	C::Vector{AbstructArray}
end

function YAlphaWaveFunc(q,rmesh,nu,ParamIndex)
	potset=LambdaAlphaPot.CalcPots(rmesh,nu,ParamIndex)
	χ=Scattering.RadWaveFunc(q,potset)
	return χ
end

function CalcCF(q,rmesh,nu,ParamIndex)
	χ=YAlphaWaveFunc(q,rmesh,nu,ParamIndex)

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