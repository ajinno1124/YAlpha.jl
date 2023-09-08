module SkyrmeParams
    using CSV, DataFrames

    function read_SkyrmeParam()
		input_file="LambdaParameters.dat"
		if isfile(input_file)==false
			println("Put LambdaParameters.dat in the run file directory.")
			#println("Using sample file in share.")
			#df_Lambda=DataFrame(CSV.File("../../share/LambdaParameters.dat", delim='\t', comment="#"))
			exit(1)
		else
        	df_Lambda = DataFrame(CSV.File(input_file, delim='\t', comment="#"))
		end
		return df_Lambda
    end

	#const df_Lambda=read_SkyrmeParam()
	df_Lambda=read_SkyrmeParam()

    function getaL(ParamIndex::Int)
        aL=zeros(Float64,5)
        args=["a1","a2","a3","a4","a5"]
        for i=eachindex(args)
            aL[i]=df_Lambda[ParamIndex,args[i]]
        end

        return aL
    end

    export read_SkyrmeParam, getaL, df_Lambda

end