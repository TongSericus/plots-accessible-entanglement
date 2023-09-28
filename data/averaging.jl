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
    filelist::Vector{Vector{String}}, path::String, filename::String
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
    filelist::Vector{String}, path::String, filename::String;
    param_name::String = "Null", param_list::Vector{Float64} = [0.0]
)
    L = length(filelist)
    S₂_avg = zeros(Float64, L)
    S₂_err = zeros(Float64, L)
    for (n, f) in enumerate(filelist)
        data = load("$(path)/$(f)")
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

# Read and process probability data to obtain Pn, Pm, and Pmn
function merge_prob_data(
    filenames::Vector{String}, path::String, filename::String;
    p::String = "Pn", isConj::Bool = true
)
    # raw data is saved in multiple files
    data = load.("$(path)" .* filenames)
    # merge raw data
    tmp = [data[i][p*"_up"] for i in 1:length(filenames)]
    Pₙ_up = hcat(tmp...)
    isConj ? Pₙ_dn = conj.(Pₙ_up) : 
        (
             tmp = [data[i][p*"_dn"] for i in 1:length(filenames)];
             Pₙ_dn = hcat(tmp...)
        )

    LA = size(Pₙ_up,1)
    L = size(Pₙ_up,2)
    Pₙ_list = zeros(ComplexF64, 2*LA-1, L)
    Pₘ_list = zeros(ComplexF64, 2*LA-1, L)
    Pₘₙ_list = zeros(ComplexF64, LA, LA, L)

    for n in axes(Pₙ_up,2)
        @views kron!(Pₘₙ_list[:, :, n], Pₙ_up[:, n], Pₙ_dn[:, n])
        @views sum_anti_diag!(Pₙ_list[:, n], Pₘₙ_list[:, :, n])
        @views sum_diag!(Pₘ_list[:, n], Pₘₙ_list[:, :, n])
    end

    p == "Pn" ? str = ["Pmn", "Pn", "Pm"] : str = ["Pmn2", "Pn2", "Pm2"]
    # save the merged data
    jldopen("$(filename)", "w") do file
        write(file, str[1], Pₘₙ_list)
        write(file, str[2], Pₙ_list)
        write(file, str[3], Pₘ_list)
    end
end

# Average the merged probability data
function average_prob_data(
    filelist::Vector{String}, LA::Int, path::String, filename::String;
    p::String = "Pn",
    param_name::String = "Null", param_list::Vector{Float64} = [0.0]
)
    L = length(filelist)
    Pₘₙ_avg = zeros(ComplexF64, LA+1, LA+1, L)
    Pₘₙ_err = zeros(ComplexF64, LA+1, LA+1, L)
    Pₙ_avg = zeros(ComplexF64, 2*LA+1, L)
    Pₙ_err = zeros(ComplexF64, 2*LA+1, L)
    Pₘ_avg = zeros(ComplexF64, 2*LA+1, L)
    Pₘ_err = zeros(ComplexF64, 2*LA+1, L)

    p == "Pn" ? str = ["Pmn", "Pn", "Pm"] : str = ["Pmn2", "Pn2", "Pm2"]

    for (n, f) in enumerate(filelist)
        data = load("$(path)/$(f)")
        # raw data is assumed to be approximately uncorrelated
        @views Pₘₙ_avg[:, :, n] .= mean(data[str[1]], dims=3)
        @views Pₘₙ_err[:, :, n] .= std(data[str[1]], dims=3) / sqrt(size(data[str[1]],3))

        @views Pₙ_avg[:, n] .= mean(data[str[2]], dims=2)
        @views Pₙ_err[:, n] .= std(data[str[2]], dims=2) / sqrt(size(data[str[2]],2))

        @views Pₘ_avg[:, n] .= mean(data[str[3]], dims=2)
        @views Pₘ_err[:, n] .= std(data[str[3]], dims=2) / sqrt(size(data[str[3]],2))
    end

    # save the processed data
    jldopen("./processed_data/$(filename)", "w") do file
        write(file, "$(param_name)", param_list)
        write(file, str[1] * "_avg", Pₘₙ_avg)
        write(file, str[1] * "_err", Pₘₙ_err)
        write(file, str[2] * "_avg", Pₙ_avg)
        write(file, str[2] * "_err", Pₙ_err)
        write(file, str[3] * "_avg", Pₘ_avg)
        write(file, str[3] * "_err", Pₘ_err)
    end
end

### Estimate with Jack's Knife ###
# estimate Shannon entropy
function estimate_Hα(Pn::AbstractArray{T}; α::Float64 = 0.5) where T
    # find a rough scope of meaningful distributions
    Pn_avg = reshape(mean(Pn, dims=2), length(axes(Pn,1)))
    idx = findall(x -> x > 1e-6, Pn_avg)

    Pn_list = zeros(T, length(idx), length(axes(Pn,2)))
    Hα_jack = zeros(T, length(axes(Pn,2)))
    for (i,n) in enumerate(idx)
        Pn_list[i, :] = JackknifeObservable(@view Pn[n, :])
    end

    for i in axes(Pn,2)
        idx = findall(x -> x > 1e-6, @view Pn_list[:, i])
        Hα_jack[i] = Hα(α, @view Pn_list[idx, i])
    end

    Hα_avg = mean(Hα_jack)
    Hα_err = std(Hα_jack) * sqrt(length(Hα_jack))

    return Hα_avg, Hα_err
end