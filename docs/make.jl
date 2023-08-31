using YAlpha
using Documenter

DocMeta.setdocmeta!(YAlpha, :DocTestSetup, :(using YAlpha); recursive=true)

makedocs(;
    modules=[YAlpha],
    authors="ajinno1124 <jinno.asanosuke.36w@st.kyoto-u.ac.jp> and contributors",
    repo="https://github.com/ajinno1124/YAlpha.jl/blob/{commit}{path}#{line}",
    sitename="YAlpha.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ajinno1124.github.io/YAlpha.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ajinno1124/YAlpha.jl",
    devbranch="main",
)
