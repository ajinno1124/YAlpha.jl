module Output

using Printf, CSV, DataFrames

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
include("./tune_a3Lam.jl")
import .tune_a3Lam

function PrintHeader(io,rmesh,nu,ParamIndex,withmom)
	println(io,"# nu = ", nu)
	println(io,"# rmesh =", rmesh)
	println(io,"# Parameter Index = ", ParamIndex)
	println(io,"# withmom = ", withmom)
end

function PrintHeader(io,rmesh,nu,ParamIndex)
	println(io,"# nu = ", nu)
	println(io,"# rmesh =", rmesh)
	println(io,"# Parameter Index = ", ParamIndex)
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


function Output_a3opt(E_ans,rmesh,nu,ParamIndex)    
    io1=open("a3.dat","w")
    PrintHeader(io1,rmesh,nu,ParamIndex)
	println(io1,"# E_ans = ",E_ans)

	println(io1,"ParameterName	a3	BE(5HeLam)(MeV)	BE-E_ans(MeV)")

    for i=eachindex(ParamIndex)
		a3_opt=tune_a3Lam.Optimize_a3(E_ans,rmesh,nu,ParamIndex[i])
		aL=SkyrmeParams.getaL(ParamIndex[i])
		BE=tune_a3Lam.Calc_BE_a3(a3_opt,rmesh,nu,aL)

		@printf(io1,"%s\t",df_Lambda[ParamIndex[i],"ParameterName"])
		@printf(io1,"%1.5f\t",a3_opt)
		@printf(io1,"%1.5f\t",BE)
		@printf(io1,"%1.5f\n",BE-E_ans)
    end
	
	close(io1)

end

function Replace_a3(E_ans,rmesh,nu,ParamIndex)
	df_a3=DataFrame(CSV.File("a3.dat", delim='\t', comment="#"))
	println(df_a3)

	io1=open("LambdaParameters_tunea3.dat","w")
	PrintHeader(io1,rmesh,nu,ParamIndex)
	println(io1,"# E_ans = ",E_ans)

	println(io1,"ParameterName	a1	a2	a3	a4	a5	BE(5He_Lam)	BE-E_ans(MeV)")

	for i=eachindex(ParamIndex)
		@printf(io1,"%s\t",df_Lambda[ParamIndex[i],"ParameterName"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"a1"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"a2"])
		@printf(io1,"%1.5f\t",df_a3[ParamIndex[i],"a3"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"a4"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"a5"])
		@printf(io1,"%1.5f\t",df_a3[ParamIndex[i],"BE(5HeLam)(MeV)"])
		@printf(io1,"%1.5f\n",df_a3[ParamIndex[i],"BE-E_ans(MeV)"])
	end
end

export Output_BoundState, Output_Potential, Output_PhaseShift, Output_CF, Output_a3opt, Replace_a3

end