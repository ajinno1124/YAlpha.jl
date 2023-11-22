module LamAlphaBoundState

include("./constants.jl")
include("./LamAlphaPot.jl")
import .LamAlphaPot
include("./Scattering.jl")
using .Scattering
include("./MyLib.jl")
import .MyLib

function BoundCond(E,mass,f_vec,rmesh)
    @assert E<=0
    h=rmesh[2]-rmesh[1]
	N_rmesh=length(rmesh)
	l=0

    Rin=zeros(Float64,3)
    Rin[1]=0.5*h # r=rmesh[1]=0.5*h
	Rin[2]=Numerov6((-1)^(l+1)*Rin[1],Rin[1], view(f_vec,vcat(1,1:2)),h)
    Rin[3]=Numerov6(view(Rin,1:2), view(f_vec,1:3),h)

    Rout=zeros(Float64,3)
    Rout[3]=exp(-(-2*mass/ħc^2*E)^(0.5) * rmesh[N_rmesh])
    Rout[2]=exp(-(-2*mass/ħc^2*E)^(0.5) * rmesh[N_rmesh-1])
    Rout[1]=exp(-(-2*mass/ħc^2*E)^(0.5) * rmesh[N_rmesh-2])

    return Rin,Rout
end

function Calc_fvec(E,μ,rmesh,PS)
	f_vec=zeros(Float64,length(rmesh))
	@. f_vec += -0.25*PS.dh2_2μeff[:]^2/PS.h2_2μeff[:] + 0.5*PS.ddh2_2μeff[:]
	@. f_vec += PS.U[:] + PS.dh2_2μeff[:]/rmesh[:]
	@. f_vec -= E
	@. f_vec /= -PS.h2_2μeff

	return f_vec
end

function WronskyEuler(E,μ,rmesh,PS)
    h=rmesh[2]-rmesh[1]
	N_rmesh=length(rmesh)
	Nmatch=floor(Int,length(rmesh)/2)

	f_vec=Calc_fvec(E,μ,rmesh,PS)

    Rin=zeros(Float64,5)
    Rout=zeros(Float64,5)
    Rin[3:5],Rout[1:3]=BoundCond(E,μ,f_vec,rmesh)


    for i in 3:Nmatch+1
		for n in 1:4
            Rin[n]=Rin[n+1]
		end
        ψvec=view(Rin,3:4)
        fvec=view(f_vec,i-1:i+1)
        Rin[5]=Numerov6(ψvec,fvec,h)
    end

    for i in N_rmesh-2:-1:Nmatch-1
		for n in 5:-1:2
            Rout[n]=Rout[n-1]
		end
        ψvec=view(Rout,3:-1:2)
        fvec=view(f_vec,i+1:-1:i-1)
        Rout[1]=Numerov6(ψvec,fvec,-h)
    end

    dRin=MyLib.diff1st5pt(h,Rin)
    dRout=MyLib.diff1st5pt(h,Rout)

    return Rin[3]*dRout-Rout[3]*dRin
end

export WronskyEuler

function FindE(μ,rmesh,PS)
	Erange=-10.01:1.0:-0.01
	args=[μ,rmesh,PS]
    for i in 1:(length(Erange)-1)
        Eans=MyLib.MyBisect(Erange[i],Erange[i+1],WronskyEuler,args,rtol=1e-6)
		if isnan(Eans)==true
            continue
        else
            return Eans
        end
	end
	return NaN
end

function NormFactor(rmesh,ψ)
    ans=MyLib.IntTrap(rmesh,@. ψ[:]^2)
	ans+=0.5*ψ[1]^2 #add between r=0~rmesh[1]
    ans=sqrt(ans)
    return ans
end

function ShootWaveFunc(E,μ,rmesh,PS)
    h=rmesh[2]-rmesh[1]
	Nmatch=floor(Int,length(rmesh)/2)

	f_vec=Calc_fvec(E,μ,rmesh,PS)

	N_rmesh=length(rmesh)
    R=zeros(Float64,N_rmesh)
    R[1:3],R[N_rmesh-2:N_rmesh]=BoundCond(E,μ,f_vec,rmesh)

    for i in 3:Nmatch-1
        ψvec=view(R,i-1:i)
        fvec=view(f_vec,i-1:i+1)
        R[i+1]=Numerov6(ψvec,fvec,h)
    end

	for i in 1:Nmatch
		R[i]/=R[Nmatch]
	end

    for i in N_rmesh-2:-1:Nmatch+1
        ψvec=view(R,i+1:-1:i)
        fvec=view(f_vec,i+1:-1:i-1)
        R[i-1]=Numerov6(ψvec,fvec,-h)
    end

	for i in Nmatch+1:N_rmesh
		R[i]/=R[Nmatch]
	end
	R[Nmatch]/=R[Nmatch]

    @. R[:]*=(PS.h2_2μeff[:])^(-0.5)

    Norm=NormFactor(rmesh,R)
    R*=sign(R[2]-R[1])/Norm

    return R
end

function BoundWaveFunc(μ,rmesh,PS::PotSet)
	h=rmesh[2]-rmesh[1]
	E=FindE(μ,rmesh,PS)
	if isnan(E)==true
		ψ=fill(NaN,length(rmesh))
        println("No bound state is found. NaN is returned.")
	else
		ψ=ShootWaveFunc(E,μ,rmesh,PS)
	end
	return E,ψ
end

function Calc_LamAlphaBoundState(rmesh,nu,ParamIndex,df_Lambda; withmom=true, Gauss=false)
	LamAPot=LamAlphaPot.CalcPotentials(rmesh,nu,ParamIndex,df_Lambda,Gauss=Gauss)
	μ=mΛMeV*mαMeV/(mΛMeV + mαMeV)
	if withmom==false
		LamAPot.h2_2μeff=ħc^2/(2*μ)*ones(Float64,length(rmesh))
		LamAPot.dh2_2μeff=zeros(Float64,length(rmesh))
		LamAPot.ddh2_2μeff=zeros(Float64,length(rmesh))
	end
	PS=Scattering.PotSet(LamAPot.h2_2μeff, LamAPot.dh2_2μeff, LamAPot.ddh2_2μeff, LamAPot.U_local)
	return BoundWaveFunc(μ,rmesh,PS)
end

function DistributionSize(rmesh::AbstractArray,ψ::AbstractArray)
    r=MyLib.IntTrap(rmesh,(@. rmesh[:]^2*ψ[:]^2))
    r+=0.5*rmesh[1]^2*ψ[1]^2
    r=sqrt(r)
    return r
end

export Calc_LamAlphaBoundState

end