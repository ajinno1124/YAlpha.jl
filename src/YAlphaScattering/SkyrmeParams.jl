module SkyrmeParams
    using CSV, DataFrames

    function read_SkyrmeParam(input_file)
        df_Lambda = DataFrame(CSV.File(input_file, delim='\t', comment="#"))
        return df_Lambda
    end

    function getaL(df_Lambda::DataFrame,ParameterIndex::Int)
        params=5
        aL=zeros(Float64,params)
        args=["a1","a2","a3","a4","a5"]
        for i=eachindex(args)
            aL[i]=df_Lambda[ParameterIndex,args[i]]
        end

        return aL
    end

    export read_SkyrmeParam, getaL

end