using YAlpha, DataFrames

function run_Bound()
	h=0.1
	N_rmesh=200
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
	nu=0.2
	qcmesh=10:200
	R=[1.0,3.0,5.0]

	ParamIndex=1:nrow(df_Lambda)
	println(df_Lambda)

	#Output_BoundState(rmesh,nu,ParamIndex,withmom=true)
	#Output_Potential(rmesh,nu,ParamIndex)
	#Output_PhaseShift(qcmesh,rmesh,nu,ParamIndex)
	Output_CF(qcmesh,rmesh,nu,ParamIndex,R)
end

run_Bound()