module tune_nu

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


function Calc_BE_nu(nu::Float64,rmesh,aL,γ)
    LamAPot=LamAlphaPot.CalcPotentials(rmesh,nu,aL,γ)
	μ=mΛMeV*mαMeV/(mΛMeV + mαMeV)
	PS=Scattering.PotSet(LamAPot.h2_2μeff, LamAPot.dh2_2μeff, LamAPot.ddh2_2μeff, LamAPot.U_local)
    E=LamAlphaBoundState.FindE(μ,rmesh,PS)
    return E
end

function dev_BE(nu,E_ans,rmesh,aL,γ)
    E_cal=Calc_BE_nu(nu,rmesh,aL,γ)
    return E_cal-E_ans
end


function Optimize_nu(E_ans,rmesh,ParamIndex::Int,df_Lambda)
    nu=0.18:0.001:0.4
    aL,γ=SkyrmeParams.getaL_gamma(df_Lambda,ParamIndex)
    args=(E_ans,rmesh,aL,γ)

    for i in 1:length(nu)-1
        nu_opt=MyLib.MyBisect(nu[i], nu[i+1], dev_BE, args)
        if isnan(nu_opt)==true
            continue
        else
            return nu_opt
        end
    end

    println("Optimized nu is not found.")
    return NaN
end

export Calc_BE_nu, dev_BE, Optimize_nu

end