module tune_a3Lam

include("./constants.jl")
include("./LamAlphaPot.jl")
import .LamAlphaPot
include("./SkyrmeParams.jl")
import .SkyrmeParams
include("./Scattering.jl")
import .Scattering
include("./LamAlphaBoundState.jl")
import .LamAlphaBoundState
include("./MyLib.jl")
import .MyLib


function Calc_BE_a3(a3::Float64,rmesh,nu,aL,γ)
    aL[3]=a3
    LamAPot=LamAlphaPot.CalcPotentials(rmesh,nu,aL,γ)
	μ=mΛMeV*mαMeV/(mΛMeV + mαMeV)
	PS=Scattering.PotSet(LamAPot.h2_2μeff, LamAPot.dh2_2μeff, LamAPot.ddh2_2μeff, LamAPot.U_local)
    E=LamAlphaBoundState.FindE(μ,rmesh,PS)
    return E
end

function dev_BE(a3,E_ans,rmesh,nu,aL,γ)
    E_cal=Calc_BE_a3(a3,rmesh,nu,aL,γ)
    return E_cal-E_ans
end


function Optimize_a3(E_ans,rmesh,nu,ParamIndex::Int,df_Lambda)
    a3=-100.0:20.0:100.0
    aL,γ=SkyrmeParams.getaL_gamma(df_Lambda,ParamIndex)
    args=(E_ans,rmesh,nu,aL,γ)

    for i in 1:length(a3)-1
        a3_opt=MyLib.MyBisect(a3[i], a3[i+1], dev_BE, args)
        if isnan(a3_opt)==true
            continue
        else
            return a3_opt
        end
    end

    println("No optimized a3 is found.")
    return NaN
end

export Calc_BE_a3, dev_BE, Optimize_a3

end