module SkyrmeParams

using CSV, DataFrames

function read_SkyrmeParam(input_file)
    df_Lambda=DataFrame[]

    if isfile(input_file)==false
        print("input file should be a valid file.")
    else
        df_Lambda = DataFrame(CSV.File(input_file, delim='\t', comment="#"))
    end

    return df_Lambda
end

function getaL_gamma(df_Lambda,ParamIndex::Int)
    aL=zeros(Float64,5)
    γ=zeros(Float64,2)
    args1=["a1","a2","a3","a4","a5"]
    args2=["gamma1","gamma2"]
    for i=eachindex(args1)
        aL[i]=df_Lambda[ParamIndex,args1[i]]
    end
    for i=eachindex(args2)
        γ[i]=df_Lambda[ParamIndex,args2[i]]
    end

    return aL,γ

end

function get_GaussVa(df_Lambda,ParamIndex::Int)
    num=3
    Vi=zeros(Float64,num)
    ai=zeros(Float64,num)
    args1=["V1","V2","V3"]
    args2=["a1","a2","a3"]
    for i=eachindex(args1)
        Vi[i]=df_Lambda[ParamIndex,args1[i]]
    end
    for i=eachindex(args2)
        ai[i]=df_Lambda[ParamIndex,args2[i]]
    end

    return Vi,ai

end

export getaL_gamma, read_SkyrmeParam, get_GaussVa


end