# Script for averaging entanglement entropy data

include("averaging.jl")

path = "./Pq2/Lx8_Ly8_LA8mul2/N32/"
U_list = [-1.0, -2.0, -3.0, -4.0, -6.0, -8.0]

# merge raw data
for U in U_list
    filenames = filter(x -> match(Regex("GS_Pn2_U$(U)_N32_Lx8_Ly8_LA16_beta18.0_L180_*"), x) !== nothing, readdir(path))
    merge_prob_data(
        filenames, path, path * "merged_data/GS_Pn2_U$(U)_N32_Lx8_Ly8_LA16_beta18.0_L180.jld",
        p = "Pn2", isConj = true
    )
end

# average merged data
filelist = ["GS_Pn2_U$(U)_N32_Lx8_Ly8_LA16_beta18.0_L180.jld" for U in U_list]
average_prob_data(
    filelist, 16, path * "merged_data", "GS_Pn2_N32_Lx8_Ly8_LA16_beta18.0_L180.jld",
    p = "Pn2", param_name = "U", param_list = U_list
)
