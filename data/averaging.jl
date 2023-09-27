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
function process_num_data(filenames::Vector{String}, path::String, filename::String)
    data = load.("$(path)/$(num_files[1])")
    
    p_num_list = data["p_num"]
    for i in 2 : length(num_files)
        data = JLD.load("$(path)/$(num_files[i])")

        p_num_list = vcat(p_num_list, data["p_num"])
    end

    jldopen("./processed_data/$(filename)", "w") do file
        write(file, "p_num", p_num_list)
    end

    return nothing
end

function process_denom_data(denom_files::Vector{String}, path::String, filename::String)

    data = JLD.load("$(path)/$(denom_files[1])")
    
    p_denom_list = data["p_denom"]
    Pn2_up_list = data["Pn2_up"]
    Pn2_dn_list = data["Pn2_dn"]

    for i in 2 : length(denom_files)
        data = JLD.load("$(path)/$(denom_files[i])")

        p_denom_list = vcat(p_denom_list, data["p_denom"])
        Pn2_up_list = hcat(Pn2_up_list, data["Pn2_up"])
        Pn2_dn_list = hcat(Pn2_dn_list, data["Pn2_dn"])
    end

    LA = length(Pn2_up_list[:, 1])
    Pmn_tmp = zeros(eltype(Pn2_up_list), LA, LA)
    nsamples = size(Pn2_up_list)[2]
    Pn_list = zeros(eltype(Pn2_up_list), 2*LA-1, nsamples)
    Pm_list = zeros(eltype(Pn2_up_list), 2*LA-1, nsamples)

    for i in 1 : nsamples
        @views kron!(Pmn_tmp, Pn2_up_list[:, i], Pn2_dn_list[:, i])
        @views sum_anti_diag!(Pn_list[:, i], Pmn_tmp)
        @views sum_diag!(Pm_list[:, i], Pmn_tmp)
    end
    
    jldopen("$(filename)", "w") do file
        write(file, "p_denom", p_denom_list)
        write(file, "Pn", real(Pn_list))
        write(file, "Pm", real(Pm_list))
    end

    return nothing
end

function process_Pn_data(Pn_files::Vector{String}, path::String, filename::String)

    data = JLD.load("$(path)/$(Pn_files[1])")
    
    Pn_up_list = data["Pn_up"]
    Pn_dn_list = data["Pn_dn"]

    for i in 2 : length(Pn_files)
        data = JLD.load("$(path)/$(Pn_files[i])")

        Pn_up_list = hcat(Pn_up_list, data["Pn_up"])
        Pn_dn_list = hcat(Pn_dn_list, data["Pn_dn"])
    end

    LA = length(Pn_up_list[:, 1])
    Pmn_tmp = zeros(eltype(Pn_up_list), LA, LA)
    nsamples = size(Pn_up_list)[2]
    Pn_list = zeros(eltype(Pn_up_list), 2*LA-1, nsamples)
    Pm_list = zeros(eltype(Pn_up_list), 2*LA-1, nsamples)

    for i in 1 : nsamples
        @views kron!(Pmn_tmp, Pn_up_list[:, i], Pn_dn_list[:, i])
        @views sum_anti_diag!(Pn_list[:, i], Pmn_tmp)
        @views sum_diag!(Pm_list[:, i], Pmn_tmp)
    end
    
    jldopen("$(filename)", "w") do file
        write(file, "Pn", real(Pn_list))
        write(file, "Pm", real(Pm_list))
    end

    return nothing
end

### Estimate with Jack's Knife ###
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

# Renyi-2 entanglement entropy
function estimate_S2(p_num::AbstractVector, p_denom::AbstractVector)
    # numerator
    p_num_jack = JackknifeObservable(p_num)

    J_num_avg = mean(p_num_jack)
    J_num_err = std(p_num_jack) * sqrt(length(p_num_jack))
    J_num = measurement(J_num_avg, J_num_err)

    # denominator
    p_denom_jack = JackknifeObservable(p_denom)

    J_denom_avg = mean(p_denom_jack)
    J_denom_err = std(p_denom_jack) * sqrt(length(p_denom_jack))
    J_denom = measurement(J_denom_avg, J_denom_err)

    return -log(J_num / J_denom)
end