using JLD, Measurements

include("./averaging_utils.jl")

### Estimators ###
# Renyi-generalized Shannon entropy
Hα(α::E, Pnα::AbstractArray{T}) where {T,E} = begin
    # regular Shannon entropy -∑p⋅lnp
    α == 1 && return -sum(Pnα .* log.(Pnα))
    α == 1.0 && return -sum(Pnα .* log.(Pnα))

    return log(sum(Pnα .^ α)) / (1 - α)
end

# Renyi-generalized symmetry-resolved entanglement
ΔSα(Pn::AbstractArray{T}, Pnα::AbstractArray, α::Int) where T = log(Pnα / Pn^α)

### Process raw data and averaging ###
# Read and process entanglement entropy data that is stored as det(g_A)
function merge_EE_data(
    filelist::Vector{Vector{String}}, path::String, filename::String;
)   
    detgA_list = Vector{Float64}[]
    for filenames in filelist
        # raw data is saved in multiple files
        data = load.("$(path)" .* filenames)
        # merge raw data
        tmp = [data[i]["detgA"] for i in 1:length(filenames)]
        push!(detgA_list, vcat(tmp...))
    end

    # save the merged data
    jldopen("$(filename)", "w") do file
        write(file, "detgA", detgA_list)
    end
end

# Average the merged EE data
function average_EE_data(
    filelists::Vector{Vector{String}}, path::String, filename::String;
    param_name::String = "Null", param_list::Vector{Float64} = [0.0]
)
    L = length(filelists)
    S₂_avg = zeros(Float64, L)
    S₂_err = zeros(Float64, L)
    for (n, filenames) in enumerate(filelists)
        data = load.("$(path)/$(filenames)")
        # raw data is assumed to be approximately uncorrelated
        tmp_avg = mean.(data["detgA"])
        tmp_err = std.(data["detgA"]) ./ length.(data["detgA"])

        # error propagation
        expmS₂ = measurement.(tmp_avg, tmp_err)
        S₂ = -log(prod(expmS₂))
        S₂_avg[n] = S₂.val
        S₂_err[n] = S₂.err
    end

    # save the processed data
    jldopen("./processed_data/$(filename)", "w") do file
        write(file, "$(param_name)", param_list)
        write(file, "S2_avg", S₂_avg)
        write(file, "S2_err", S₂_err)
    end
end
