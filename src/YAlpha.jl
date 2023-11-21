module YAlpha

    include("./YAlphaScattering/constants.jl")
    include("./YAlphaScattering/SkyrmeParams.jl")
    using .SkyrmeParams
    export read_SkyrmeParam,getaL_gamma,df_Lambda,ħc

    include("./YAlphaScattering/LamAlphaPot.jl")
    using .LamAlphaPot
    export CalcPotentials, GaussianPot

	#for debugging
	include("./YAlphaScattering/Scattering.jl")
	using .Scattering
	export RadWaveFunc, PotSet, PhaseShift

	include("./YAlphaScattering/CorrelationFunc.jl")
	using .CorrelationFunc
	export LamAlphaWaveFunc, CorrelationFunction, LLformula_S1

	include("./YAlphaScattering/LamAlphaBoundState.jl")
	using .LamAlphaBoundState
	export Calc_LamAlphaBoundState, WronskyEuler

	include("./YAlphaScattering/Output.jl")
	using .Output
	export Output_BoundState, Output_Potential, Output_PhaseShift, Output_CF, Output_a3opt, Replace_a3, Output_LL
	export Output_nuopt

	include("./Fitting/FittingSkyrme.jl")
	using .FittingSkyrme
	export Fit_Ulocal, OutputFitResult

end
