using YAlpha

function run_Bound()
	h=0.1
	N_rmesh=500
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
	nu=0.2
	ParamIndex=[9,10,11]
	println(df_Lambda)

	Output_BoundState(rmesh,nu,ParamIndex)
end

run_Bound()