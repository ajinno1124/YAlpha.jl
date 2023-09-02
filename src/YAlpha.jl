module YAlpha

    include("./YAlphaScattering/constants.jl")
    include("./YAlphaScattering/SkyrmeParams.jl")
    using .SkyrmeParams
    export read_SkyrmeParam,getaL,df_Lambda,ħc

    include("./YAlphaScattering/LamAlphaPot.jl")
    using .LamAlphaPot
    export CalcPotentials

	include("./YAlphaScattering/CorrelationFunc.jl")
	using .CoorelationFunc
	export LamAlphaWaveFunc

end
