module SkyrmeParams

using CSV, DataFrames

function read_SkyrmeParam(input_file)
    #input_file="LambdaParameters.dat"
    df_Lambda=DataFrame[]

    if isfile(input_file)==false
        print("input file should be a valid file.")
    else
        df_Lambda = DataFrame(CSV.File(input_file, delim='\t', comment="#"))
    end

    return df_Lambda
end

#df_Lambda=read_SkyrmeParam(input_file)
#df_Lambda=read_SkyrmeParam(ARGS)

function getaL(df_Lambda,ParamIndex::Int)
    aL=zeros(Float64,5)
    args=["a1","a2","a3","a4","a5"]
    for i=eachindex(args)
        aL[i]=df_Lambda[ParamIndex,args[i]]
    end

    return aL

end

export getaL, read_SkyrmeParam


end