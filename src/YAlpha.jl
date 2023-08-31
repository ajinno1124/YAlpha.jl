module YAlpha

    include("YAlphaScattering/constants.jl")
    include("YAlphaScattering/SkyrmeParams.jl")
    using .SkyrmeParams
    export read_SkyrmeParam,getaL, ħc

    include("YAlphaScattering/LamAlphaPot.jl")
    using .LamAlphaPot
    export CalcPotentials

end
