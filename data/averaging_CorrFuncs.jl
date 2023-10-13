# Script for averaging entanglement entropy data

include("averaging.jl")

path = "./CorrFuncs/N32/"
U_list = collect(-1.0:-1.0:-8.0)

# merge raw data
for U in U_list
    filenames = filter(x -> match(Regex("GS_CorrFuncs_U$(U)_N32_Lx8_Ly8_beta18.0_*"), x) !== nothing, readdir(path))
    merge_corrfuncs_data(
        filenames, path, path * "merged_data/GS_PairCorr_N32_U$(U)_Lx8_Ly8_beta18.0.jld";
        param = "Ps"
    )
end

# average merged data
filelist = ["GS_PairCorr_N32_U$(U)_Lx8_Ly8_beta18.0.jld" for U in U_list]
average_corrfuncs_data(
    filelist, path * "merged_data/", "GS_PairCorr_Lx8_Ly8_N32_beta18.0.jld";
    param = "Ps", q = (0,0)
)
