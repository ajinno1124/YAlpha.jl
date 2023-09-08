module Output

using CSV, DataFrames
include("./LamAlphaBoundState.jl")
using .LamAlphaBoundState
include("./SkyrmeParams.jl")
using .SkyrmeParams

function PrintHeader(io,rmesh,nu,ParamIndex,withmom)
	println(io,"# nu = ", nu)
	println(io,"# Parameter Index = ", ParamIndex)
	println(io,"# withmom = ", withmom)
end

function Output_BoundState(rmesh,nu,ParamIndex,withmom=true)
	file_path="data/BoundState"
	rm(file_path,force=true,recursive=true)
	mkpath(file_path)
	io1=open("$(file_path)/BindingEnergy.dat","w")
	PrintHeader(io1,rmesh,nu,ParamIndex,withmom)
	println(io1,"ParameterName B.E.(MeV)")

	for i=eachindex(ParamIndex)
		io2=open("$(file_path)/$(df_Lambda[ParamIndex[i],"ParameterName"]).dat","w")
		PrintHeader(io2,rmesh,nu,ParamIndex[i],withmom)
		println(io2,"r(fm) u")

		E,ψ=Calc_LamAlphaBoundState(rmesh,nu,ParamIndex[i],withmom=true)
		println(io1,df_Lambda[ParamIndex[i],"ParameterName"], " ", E)
		for j=eachindex(rmesh)
			println(io2,rmesh[j], " ", ψ[j])
		end
		close(io2)
	end
	close(io1)

end

export Output_BoundState

end