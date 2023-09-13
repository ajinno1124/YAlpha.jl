using YAlpha, DataFrames

function run_Bound()
	h=0.1
	N_rmesh=200
	rmesh=0.5*h:h:h*(N_rmesh-0.5)
	nu=0.27
	qcmesh=10:200
	R=[1.0,3.0,5.0]
	input_file="./LambdaParameters.dat"

	E_ans = -3.12
	ParamIndex=1:8
	Output_a3opt(E_ans,rmesh,nu,ParamIndex,input_file)
	Replace_a3(E_ans,rmesh,nu,ParamIndex,input_file)

	input_file2="./LambdaParameters_tunea3.dat"
	println(read_SkyrmeParam(input_file2))

	Output_BoundState(rmesh,nu,ParamIndex,input_file2,withmom=true)
	Output_Potential(rmesh,nu,ParamIndex,input_file2)
	Output_PhaseShift(qcmesh,rmesh,nu,ParamIndex,input_file2)
	Output_CF(qcmesh,rmesh,nu,ParamIndex,R,input_file2)

end

run_Bound()