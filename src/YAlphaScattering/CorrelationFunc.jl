module LamAlphaScat

include("./Scattering.jl")
include("./LambdaAlphaPot.jl")

function YAlphaWaveFunc(rmesh,nu,ParamIndex)
	PotSet=CalcPots() #LambdaAlphaPot.jl
	

end

function CoorelationFunction()
end

end