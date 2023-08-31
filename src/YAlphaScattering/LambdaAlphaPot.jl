module LambdaAlphaPot
    include("./SkyrmeParams.jl")
    import .SkyrmeParams
    include("./constants.jl")
    include("./MyLib.jl")

    mutable struct YAlphaPotentials{T <: AbstractFloat}
        rmesh::AbstractArray
        ρ::Vector{T}
        τ::Vector{T}
        Lapρ::Vector{T}
        Ulocal::Vector{T}

        h2_2μeff::Vector{T}
        dh2_2μeff::Vector{T}
        ddh2_2μeff::Vector{T}
    end

    function GetNuc(nu::Float64,A)
        return nu/(1-1/A)
    end

    function Density(r,nuc::Float64,A)
        return (@. A*(2*nuc/π)^(3/2)*exp(-2*nuc*r^2))
    end

    function KinDensity(r,ρ,nu::Float64,A)
        return (@. (4*nu^2*r^2 + 3*nu/A)*ρ)
    end

    function LapDensity(r,ρ,nuc::Float64)
        return (@. -4*nuc*(3-4*nuc*r^2)*ρ)
    end

    function U_local(ρ,τ,Lapρ,aL)
        U=zeros(Float64,length(ρ))
        @. U=aL[1]*ρ + aL[2]*τ - aL[3]*Lapρ + aL[4]*ρ^(4/3) + aL[5]*ρ^(5/3)
        return U
    end

    function meff_Lam(ρ,a2::Float64)
        return (@. mΛMeV/(1+2*mΛMeV*a2*ρ/ħc^2))
    end

    function Calc_h2_2μeff(ρ,a2::Float64)
        meffLam=meff_Lam(ρ,a2)
        return (@. meffLam*mαMeV/(meffLam + mαMeV))
    end

    function CalcPotentials(rmesh::AbstractArray,nu,aL)
        h=rmesh[2]-rmesh[1]
        @assert h/2-rmesh[1]<0.00000001

        A=4
        nuc=GetNuc(nu,A)
        ρ=Density(rmesh,nuc,A)
        τ=KinDensity(rmesh,ρ,nu,A)
        Lapρ=LapDensity(rmesh,ρ,nuc)
        Ulocal=U_local(ρ,τ,Lapρ,aL)

        h2_2μeff=Calc_h2_2μeff(ρ,aL[2])
        dh2_2μeff=MyLib.diff1st5pt(h,h2_2μeff,1)
        ddh2_2μeff=MyLib.diff2nd5pt(h,h2_2μeff,1)

        return YAlphaPotentials(rmesh,ρ,τ,Lapρ,Ulocal,h2_2μeff,dh2_2μeff,ddh2_2μeff)
    end

    export CalcPotentials

end