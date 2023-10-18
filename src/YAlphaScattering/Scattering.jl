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

export PotSet

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

export Numerov6

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

export Calc_fvec

function normfactor(st,rmesh)
	reA=Re_A(st,rmesh)
	imA=Im_A(st,rmesh)
	return 2*(reA^2+imA^2)^(0.5)

	#reB=Re_B(st,rmesh)
	#imB=Im_B(st,rmesh)
	#return (reB^2+imB^2)^(0.5)
end

#Schrodinger eq. [-∇⋅(ħ^2/2μ*)∇ + U] ψ= Eψ, E=ħ^2*q^2/(2μ)
#Calculate wavefunction
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
	st=State(qc,R)
	norm=normfactor(st,rmesh)
	st.ψ/=norm

    return st
end

export RadWaveFunc


function Re_A(st,rmesh)
	N_rmesh=length(rmesh)
	h=rmesh[2]-rmesh[1]

	val  = st.ψ[N_rmesh]*sin(st.qc/ħc*rmesh[N_rmesh-1])
	val -= st.ψ[N_rmesh-1]*sin(st.qc/ħc*rmesh[N_rmesh])
	val /= -2*sin(st.qc/ħc*h)
	return val
end


function Im_A(st,rmesh)
	N_rmesh=length(rmesh)
	h=rmesh[2]-rmesh[1]

	val  = st.ψ[N_rmesh]*cos(st.qc/ħc*rmesh[N_rmesh-1])
	val -= st.ψ[N_rmesh-1]*cos(st.qc/ħc*rmesh[N_rmesh])
	val /= 2*sin(st.qc/ħc*h)
	return val
end

function Re_B(st,rmesh)
	N_rmesh=length(rmesh)
	h=rmesh[2]-rmesh[1]

	val  = -st.ψ[N_rmesh]*sin(st.qc/ħc*rmesh[N_rmesh-1])
	val += st.ψ[N_rmesh-1]*sin(st.qc/ħc*rmesh[N_rmesh])
	val /= 2*sin(st.qc/ħc*h)
	return val
end


function Im_B(st,rmesh)
	N_rmesh=length(rmesh)
	h=rmesh[2]-rmesh[1]

	val  = -st.ψ[N_rmesh]*cos(st.qc/ħc*rmesh[N_rmesh-1])
	val -= st.ψ[N_rmesh-1]*cos(st.qc/ħc*rmesh[N_rmesh])
	val /= 2*sin(st.qc/ħc*h)
	return val
end

function PhaseShift(st,rmesh::AbstractArray)
	#N_rmesh=length(rmesh)
	#return (asin(st.ψ[N_rmesh])-st.qc/ħc*rmesh[N_rmesh])/π

	reA=Re_A(st,rmesh)
	imA=Im_A(st,rmesh)
	return atan(reA/imA)
end

function EffRangeExp(delta,qcmesh)
	q1=qcmesh[1]
	q2=qcmesh[2]
	d1=delta[1]
	d2=delta[2]

	a0   = (q1^2-q2^2)/(q2^2*q1*cot(d1) - q1^2*q2*cot(d2))
	reff = 2*(q2*cot(d2) - q1*cot(d1))/(q2^2-q1^2)

	a0   *= ħc
	reff *= ħc

	return a0,reff
end

export PhaseShift

end