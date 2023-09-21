module FittingSkyrme

using LsqFit
include("../YALphaScattering/LamAlphaPot.jl")
import .LamAlphaPot


function Calc_args(rmesh,nu)
	h=rmesh[2]-rmesh[1]
	@assert h/2-rmesh[1]<0.00000001

	A=4
	nuc=LamAlphaPot.GetNuc(nu,A)
	ρ=LamAlphaPot.Density(rmesh,nuc,A)
	τ=LamAlphaPot.KinDensity(rmesh,ρ,nu,A)
	Lapρ=LamAlphaPot.LapDensity(rmesh,ρ,nuc)
	γ=[1/3,2/3]

	return [ρ,τ,Lapρ,γ]
end

function Calc_U_local_woMom(i::Int,aL,ρ,τ,Lapρ,γ)
	U=0.0
	U+=aL[1]*ρ[i] - aL[2]*Lapρ[i] + aL[3]*ρ[i]^(1+γ[1]) + aL[4]*ρ[i]^(1+γ[2]) #without aL2
	return U
end

function ExtractDat(pot_dat::String)
	df=DataFrame(CSV.File(pot_dat, delim='\t', comment="#"))
	return df[!,1],df[!,2]
end

function Fit_Ulocal(nu,pot_dat)
	xdata,ydata=ExtractDat(pot_dat)
	args=Calc_args(xdata,nu)
	aL0=[1.0,1.0,1.0,1.0]
	func(i,aL4)=Calc_U_local_woMom(i,aL4,args...)

	fit=curve_fit(func,eachindex(xdata),ydata,aL0)
	return fit
end

function OutputFitResult(nu,pot_dat)
	fit=Fit_Ulocal(nu,pot_dat)

	io1=open("FitResult.dat","w")
	println(io1,"# nu = $(nu)")
	println(io1,"a1	a2	a3	a4	a5")
	@printf(io1,"%.5f\t",coef(fit)[1])
	@printf(io1,"%.5f\t",0.0)
	@printf(io1,"%.5f\t",coef(fit)[2])
	@printf(io1,"%.5f\t",coef(fit)[3])
	@printf(io1,"%.5f\n",coef(fit)[4])
	close(io1)
end

export Fit_Ulocal, OutputFitResult


end