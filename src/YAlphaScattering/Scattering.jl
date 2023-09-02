module Scattering

#include("./LamAlphaPot.jl")
#import .LamAlphaPot
include("./constants.jl")

mutable struct State
    qc::AbstractFloat
    ψ::Vector{AbstractFloat}
end

mutable struct PotSet
    h2_2μeff::Vector{AbstractFloat}
	dh2_2μeff::Vector{AbstractFloat}
	ddh2_2μeff::Vector{AbstractFloat}
    U::Vector{AbstractFloat}
end

#######################################
# Differential Equation AR''+ CR=ER #
# R=u, R'=v                           #
#######################################

#Diff. eq. ψ''(r)+f(r)ψ(r)=0
#Calc ψ[1]=ψ(r-h) by ψ[2]=ψ(r) and ψ[3]=ψ(r+h) using Numerov Method
function Numerov6(ψ::AbstractArray,f::AbstractArray,h::Float64)
	val=0.0
	val+=(2-5*h^2*f[2]/6)*ψ[2]
	val-=(1+h^2*f[1]/12)*ψ[1]
	val/=(1+h^2*f[3]/12)
	return val
end

function Numerov6(ψ1::Float64,ψ2::Float64,f,h)
	val=0.0
	val+=(2-5*h^2*f[2]/6)*ψ2
	val-=(1+h^2*f[1]/12)*ψ1
	val/=(1+h^2*f[3]/12)
	return val
end

function InitialCondition!(h,R,f_vec)
    R[1]=0.5*h
	# only consider s-wave
	l=0
	R[2]=Numerov6((-1)^(l+1)*R[1],R[1], view(f_vec,vcat(1,1:2)),h)
    R[3]=Numerov6(view(R,1:2), view(f_vec,1:3),h)
end

function Calc_fvec(qc,μ,rmesh,PS::PotSet)
	f_vec=zeros(Float64,length(rmesh))
	@. f_vec += -0.25*PS.dh2_2μeff[:]^2/PS.h2_2μeff[:] + 0.5*PS.ddh2_2μeff[:]
	@. f_vec += PS.U[:] + PS.dh2_2μeff[:]/rmesh[:]
	@. f_vec -= qc^2/(2*μ)
	@. f_vec /= -PS.h2_2μeff

	return f_vec
end

#Schrodinger eq. [-∇⋅(ħ^2/2μ*)∇ + U] ψ= Eψ, E=ħ^2*q^2/(2μ)
#Calculate wavefunction and PhaseShift
function RadWaveFunc(qc,μ,rmesh,PS::PotSet)
    h=rmesh[2]-rmesh[1]
	N_rmesh=length(rmesh)

	f_vec=Calc_fvec(qc,μ,rmesh,PS)

    R=zeros(Float64,N_rmesh)
	InitialCondition!(h,R,f_vec)

    for i in 3:N_rmesh-1
        R[i+1]=Numerov6(view(R,i-1:i),view(f_vec,i-1:i+1),h)
    end

	@. R[:]*=(PS.h2_2μeff[:])^(-0.5)

	#Normalization is needed.
	#stupid method... may be able to be more faster.
	norm=findmax(view(R,floor(Int,N_rmesh/2):N_rmesh))[1]
	@. R[:]/=norm

    return State(qc,R)
end

export RadWaveFunc

function PhaseShift()
end

end