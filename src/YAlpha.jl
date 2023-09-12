module YAlpha

    include("./YAlphaScattering/constants.jl")
    include("./YAlphaScattering/SkyrmeParams.jl")
    using .SkyrmeParams
    export read_SkyrmeParam,getaL,df_Lambda,ħc

    include("./YAlphaScattering/LamAlphaPot.jl")
    using .LamAlphaPot
    export CalcPotentials

	#for debugging
	include("./YAlphaScattering/Scattering.jl")
	using .Scattering
	export RadWaveFunc, PotSet, PhaseShift

	include("./YAlphaScattering/CorrelationFunc.jl")
	using .CorrelationFunc
	export LamAlphaWaveFunc, CoorelationFunction

	include("./YAlphaScattering/LamAlphaBoundState.jl")
	using .LamAlphaBoundState
	export Calc_LamAlphaBoundState, WronskyEuler

	include("./YAlphaScattering/Output.jl")
	using .Output
	export Output_BoundState, Output_Potential, Output_PhaseShift, Output_CF, Output_a3opt, Replace_a3


end
