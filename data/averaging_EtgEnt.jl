# Script for averaging entanglement entropy data

include("averaging.jl")

path = "./EtgEnt/LA_8mul2/"
U_list = vcat(-0.5, collect(-1.0:-1.0:-8.0))
lambda = collect(0.0:0.2:0.8)

# merge raw data
for U in U_list
    filelist = [filter(x -> match(Regex("GS_attrHub_Lx8_Ly8_LA16_N32_U$(U)_lambda$(i)_beta50.0_*"), x) !== nothing, readdir(path)) for i in lambda]
    merge_EE_data(filelist, path, path * "merged_data/GS_attrHub_Lx8_Ly8_LA16_N32_U$(U)_beta50.0.jld")
end

# average merged data
filelist = ["GS_attrHub_Lx8_Ly8_LA16_N32_U$(U)_beta50.0.jld" for U in U_list]
average_EE_data(
    filelist, path * "merged_data/", "GS_EtgEnt_Lx8_Ly8_LA16_N32_beta50.0.jld";
    param_name = "U", param_list = U_list
)