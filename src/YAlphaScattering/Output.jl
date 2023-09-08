module Output

using CSV, DataFrames
include("./LamAlphaBoundState.jl")
using .LamAlphaBoundState
include("./SkyrmeParams.jl")
using .SkyrmeParams
include("./LamAlphaPot.jl")
import .LamAlphaPot
include("./Scattering.jl")
import .Scattering
include("./CorrelationFunc.jl")
import .CorrelationFunc

function PrintHeader(io,rmesh,nu,ParamIndex,withmom)
	println(io,"# nu = ", nu)
	println(io,"# Parameter Index = ", ParamIndex)
	println(io,"# withmom = ", withmom)
end

function Output_BoundState(rmesh,nu,ParamIndex;withmom=true)
	file_path="data/BoundState"
	rm(file_path,force=true,recursive=true)
	mkpath(file_path)
	io1=open("$(file_path)/BindingEnergy.dat","w")
	PrintHeader(io1,rmesh,nu,ParamIndex,withmom)
	println(io1,"ParameterName	B.E.(MeV)")

	for i=eachindex(ParamIndex)
		io2=open("$(file_path)/$(df_Lambda[ParamIndex[i],"ParameterName"]).dat","w")
		PrintHeader(io2,rmesh,nu,ParamIndex[i],withmom)
		println(io2,"r(fm)	u")

		E,ψ=Calc_LamAlphaBoundState(rmesh,nu,ParamIndex[i],withmom=true)
		println(io1,df_Lambda[ParamIndex[i],"ParameterName"], "	", E)
		for j=eachindex(rmesh)
			println(io2,rmesh[j], "	", ψ[j])
		end
		close(io2)
	end
	close(io1)

end

function Output_Potential(rmesh,nu,ParamIndex)
	file_path="data/Potentials"
	rm(file_path,force=true,recursive=true)
	mkpath(file_path)

	for i=eachindex(ParamIndex)
		io1=open("$(file_path)/Potential_$(df_Lambda[i,"ParameterName"]).dat","w")
		println(io1,"# nu = $(nu)")
		println(io1,"r(fm)	U_local(MeV)	U_m(MeV)")

		PS=LamAlphaPot.CalcPotentials(rmesh,nu,ParamIndex[i])
		U_m=zeros(Float64,length(rmesh))
		@. U_m += -0.25*PS.dh2_2μeff[:]^2/PS.h2_2μeff[:] + 0.5*PS.ddh2_2μeff[:]
		@. U_m += PS.dh2_2μeff[:]/rmesh[:]

		for j=eachindex(rmesh)
			println(io1,rmesh[j], "	", PS.U_local[j], "	", U_m[j])
		end

		close(io1)
	end
end


function Output_PhaseShift(qcmesh,rmesh,nu,ParamIndex)
	file_path="data/PhaseShift"
	rm(file_path,force=true,recursive=true)
	mkpath(file_path)

	for i=eachindex(ParamIndex)
		io1=open("$(file_path)/PhaseShift_$(df_Lambda[ParamIndex[i],"ParameterName"]).dat","w")
		println(io1,"# nu = $(nu)")
		println(io1,"q(MeV/c)	delta")

		for j=eachindex(qcmesh)
			state=CorrelationFunc.LamAlphaWaveFunc(qcmesh[j],rmesh,nu,ParamIndex[i])
			delta=Scattering.PhaseShift(state,rmesh)
			println(io1,qcmesh[j], "	",delta)
		end

		close(io1)
	end

end

function Output_CF(qcmesh,rmesh,nu,ParamIndex,R)
	file_path="data/CorrelationFunction"
	rm(file_path,force=true,recursive=true)
	mkpath(file_path)

	for i=eachindex(ParamIndex)
		for r in R
			io1=open("$(file_path)/$(df_Lambda[ParamIndex[i],"ParameterName"])_R$(r).dat","w")
			println(io1,"# nu = $(nu)")
			println(io1,"q(MeV/c)	CF")

			C=CorrelationFunc.CoorelationFunction(qcmesh,rmesh,nu,ParamIndex[i],r)

			for j=eachindex(qcmesh)
				println(io1,qcmesh[j], "	", C[j])
			end

			close(io1)
		end
	end

end

export Output_BoundState, Output_Potential, Output_PhaseShift, Output_CF

end