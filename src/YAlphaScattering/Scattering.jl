module Scattering

mutable struct State
    q::AbstractFloat
    ψ::Vector{AbstractFloat}
end

mutable struct PotSet
    h2_2μeff::Vector{AbstractFloat}
    U::Vector{AbstractFloat}
end

#######################################
# Differential Equation AR''+ CR=ER #
# R=u, R'=v                           #
#######################################

#Diff. eq. ψ''(r)+f(r)ψ(r)=0
#Calc ψ[1]=ψ(r-h) by ψ[2]=ψ(r) and ψ[3]=ψ(r+h) using Numerov Method
function Numerov6(ψ,f,h)
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

function InitialCondition()
    Rin=zeros(Float64,3)
    Rin[1]=0.5*h
	Rin[2]=Numerov6((-1)^(l+1)*Rin[1],Rin[1], view(Val,vcat(1,1:2)),h)
    Rin[3]=Numerov6(view(Rin,1:2), view(Val,1:3),h)

    return Rin,Rout
end

#Schrodinger eq. [-∇⋅(ħ^2/2μ*)∇ + U] ψ= Eψ
#Calculate wavefunction and PhaseShift
function WaveFunction(q,potset::PotSet)
    
    return ψ
end

function PhaseShift()
end

end