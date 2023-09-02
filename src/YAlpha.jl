module YAlpha

    include("YAlphaScattering/constants.jl")
    include("YAlphaScattering/SkyrmeParams.jl")
    using .SkyrmeParams
    export read_SkyrmeParam,getaL,df_Lambda,ħc

    include("YAlphaScattering/LambdaAlphaPot.jl")
    using .LambdaAlphaPot
    export CalcPotentials

end
