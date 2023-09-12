using YAlpha, DataFrames

function run_tune_a3()
	Output_a3opt(E_ans,rmesh,nu,ParamIndex)
end

function run_Bound()
	h=0.1
	N_rmesh=200
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
	nu=0.27
	qcmesh=10:200
	R=[1.0,3.0,5.0]

	ParamIndex=1:nrow(df_Lambda)
	println(df_Lambda)

	E_ans = -3.12
	Output_a3opt(E_ans,rmesh,nu,ParamIndex)
	Replace_a3(E_ans,rmesh,nu,ParamIndex)

	#Output_BoundState(rmesh,nu,ParamIndex,withmom=true)
	#Output_Potential(rmesh,nu,ParamIndex)
	#Output_PhaseShift(qcmesh,rmesh,nu,ParamIndex)
	#Output_CF(qcmesh,rmesh,nu,ParamIndex,R)

end

run_Bound()