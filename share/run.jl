using YAlpha, DataFrames

function run_Bound()
	h=0.1
	N_rmesh=500
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
	nu=0.2
	ParamIndex=1:nrow(df_Lambda)
	println(df_Lambda)

	Output_BoundState(rmesh,nu,ParamIndex,withmom=false)
	Output_Potential(rmesh,nu,ParamIndex)
end

run_Bound()