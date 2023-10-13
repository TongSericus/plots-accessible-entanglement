# Script for averaging entanglement entropy data

include("averaging.jl")

path = "./CorrFuncs/N32/"
U_list = collect(-1.0:-1.0:-8.0)

# merge raw data
for U in U_list
    filenames = filter(x -> match(Regex("GS_CorrFuncs_U$(U)_N32_Lx8_Ly8_beta18.0_*"), x) !== nothing, readdir(path))
    average_corrfuncs_data(
        filenames, path, "GS_CorrFuncs_N32_Lx8_Ly8_beta18.0.jld";
        param_list = ["Ps", "ninj_updn"]
    )
end
