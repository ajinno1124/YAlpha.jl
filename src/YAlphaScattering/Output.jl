module Output

using Printf, CSV, DataFrames

using CSV, DataFrames
include("./LamAlphaBoundState.jl")
import .LamAlphaBoundState
include("./SkyrmeParams.jl")
import .SkyrmeParams
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
	@printf(io,"# rmesh = %1.3f:%1.3f:%1.3f\n", rmesh[1],rmesh[2]-rmesh[1],rmesh[length(rmesh)])
	println(io,"# Parameter Index = ", ParamIndex)
	println(io,"# withmom = ", withmom)
end

function PrintHeader(io,rmesh,nu,ParamIndex)
	println(io,"# nu = ", nu)
	@printf(io,"# rmesh = %1.3f:%1.3f:%1.3f\n", rmesh[1],rmesh[2]-rmesh[1],rmesh[length(rmesh)])
	println(io,"# Parameter Index = ", ParamIndex)
end

function Output_BoundState(rmesh,nu,ParamIndex,input_file;withmom=true, Gauss=false)
	df_Lambda=SkyrmeParams.read_SkyrmeParam(input_file)

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

		E,ψ=LamAlphaBoundState.Calc_LamAlphaBoundState(rmesh,nu,ParamIndex[i],df_Lambda,withmom=true,Gauss=Gauss)
		r=LamAlphaBoundState.DistributionSize(rmesh,ψ)
		println(io1,"# distribution size r = $r fm")
		println(io1,df_Lambda[ParamIndex[i],"ParameterName"], "	", E)
		for j=eachindex(rmesh)
			println(io2,rmesh[j], "	", ψ[j])
		end
		close(io2)
	end
	close(io1)

end

function Output_Potential(rmesh,nu,ParamIndex,input_file; Gauss=false)
	df_Lambda=SkyrmeParams.read_SkyrmeParam(input_file)

	file_path="data/Potentials"
	rm(file_path,force=true,recursive=true)
	mkpath(file_path)

	for i=eachindex(ParamIndex)
		io1=open("$(file_path)/Potential_$(df_Lambda[i,"ParameterName"]).dat","w")
		println(io1,"# nu = $(nu)")
		println(io1,"r(fm)	U_local(MeV)	U_m(MeV)	h2_2mueff(MeV)")

		PS=LamAlphaPot.CalcPotentials(rmesh,nu,ParamIndex[i],df_Lambda, Gauss=Gauss)
		U_m=zeros(Float64,length(rmesh))
		@. U_m += -0.25*PS.dh2_2μeff[:]^2/PS.h2_2μeff[:] + 0.5*PS.ddh2_2μeff[:]
		@. U_m += PS.dh2_2μeff[:]/rmesh[:]

		for j=eachindex(rmesh)
			println(io1,rmesh[j], "	", PS.U_local[j], "	", U_m[j],"	",PS.h2_2μeff[j])
		end

		close(io1)
	end
end


function Output_PhaseShift(qcmesh,rmesh,nu,ParamIndex,input_file; Gauss=false)
	df_Lambda=SkyrmeParams.read_SkyrmeParam(input_file)

	file_path="data/PhaseShift"
	rm(file_path,force=true,recursive=true)
	mkpath(file_path)

	for i=eachindex(ParamIndex)
		io1=open("$(file_path)/PhaseShift_$(df_Lambda[ParamIndex[i],"ParameterName"]).dat","w")
		println(io1,"# nu = $(nu)")
		println(io1,"q(MeV/c)	delta")

		for j=eachindex(qcmesh)
			state=CorrelationFunc.LamAlphaWaveFunc(qcmesh[j],rmesh,nu,ParamIndex[i],df_Lambda, Gauss=Gauss)
			delta=Scattering.PhaseShift(state,rmesh)
			println(io1,qcmesh[j], "	",delta)
		end

		close(io1)
	end

end

function Output_CF(qcmesh,rmesh,nu,ParamIndex,R,input_file; Gauss=false)
	df_Lambda=SkyrmeParams.read_SkyrmeParam(input_file)

	file_path="data/CorrelationFunction"
	rm(file_path,force=true,recursive=true)
	mkpath(file_path)

	for i=eachindex(ParamIndex)
		for r in R
			io1=open("$(file_path)/$(df_Lambda[ParamIndex[i],"ParameterName"])_R$(r).dat","w")
			println(io1,"# nu = $(nu)")
			println(io1,"q(MeV/c)	CF")

			C=CorrelationFunc.CoorelationFunction(qcmesh,rmesh,nu,ParamIndex[i],r,df_Lambda,Gauss=Gauss)

			for j=eachindex(qcmesh)
				println(io1,qcmesh[j], "	", C[j])
			end

			close(io1)
		end
	end

end


function Output_a3opt(E_ans,rmesh,nu,ParamIndex,input_file)
	df_Lambda=SkyrmeParams.read_SkyrmeParam(input_file)
	#println("df_Lambda for a3")
	#println(df_Lambda)

    io1=open("a3.dat","w")
    PrintHeader(io1,rmesh,nu,ParamIndex)
	println(io1,"# E_ans = ",E_ans)

	println(io1,"ParameterName	a3	BE(5HeLam)(MeV)	BE-E_ans(MeV)")

    for i=eachindex(ParamIndex)
		a3_opt=tune_a3Lam.Optimize_a3(E_ans,rmesh,nu,ParamIndex[i],df_Lambda)
		aL,γ=SkyrmeParams.getaL_gamma(df_Lambda,ParamIndex[i])
		BE=tune_a3Lam.Calc_BE_a3(a3_opt,rmesh,nu,aL,γ)

		@printf(io1,"%s\t",df_Lambda[ParamIndex[i],"ParameterName"])
		@printf(io1,"%1.5f\t",a3_opt)
		@printf(io1,"%1.5f\t",BE)
		@printf(io1,"%1.5f\n",BE-E_ans)
    end
	
	close(io1)

end

function Replace_a3(E_ans,rmesh,nu,ParamIndex,input_file)
	df_Lambda=SkyrmeParams.read_SkyrmeParam(input_file)
	df_a3=DataFrame(CSV.File("a3.dat", delim='\t', comment="#"))
	println(df_a3)

	io1=open("LambdaParameters_tunea3.dat","w")
	PrintHeader(io1,rmesh,nu,ParamIndex)
	println(io1,"# E_ans = ",E_ans)

	println(io1,"ParameterName	a1	a2	a3	a4	a5	gamma1	gamma2	BE(5He_Lam)	BE-E_ans(MeV)")

	for i=eachindex(ParamIndex)
		@printf(io1,"%s\t",df_Lambda[ParamIndex[i],"ParameterName"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"a1"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"a2"])
		@printf(io1,"%1.5f\t",df_a3[ParamIndex[i],"a3"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"a4"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"a5"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"gamma1"])
		@printf(io1,"%1.5f\t",df_Lambda[ParamIndex[i],"gamma2"])
		@printf(io1,"%1.5f\t",df_a3[ParamIndex[i],"BE(5HeLam)(MeV)"])
		@printf(io1,"%1.5f\n",df_a3[ParamIndex[i],"BE-E_ans(MeV)"])
	end

	close(io1)
end

export Output_BoundState, Output_Potential, Output_PhaseShift, Output_CF, Output_a3opt, Replace_a3

end