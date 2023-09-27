### Utils for Data Processing ###

using Statistics

"""
    JackknifeObservable(x)

    Generate jackknife averages of the sample x
"""
function JackknifeObservable(x::AbstractVector{T}; J = zeros(T, length(x))) where T
    length(x) > 1 || throw(ArgumentError("The sample must have size > 1"))

    x_sum = sum(x)
    l = length(x)

    for i in eachindex(x)
        J[i] = (x_sum - x[i]) / (l - 1)
    end

    return J
end

"""
    kron!(C, a, b)

    Compute C = a * transpose(b) in-place
"""
function Base.kron!(C::AbstractMatrix, a::AbstractVector, b::AbstractVector)
    for i in eachindex(a)
        for j in eachindex(a)
            @inbounds C[j, i] = a[j] * b[i]
        end
    end
end

"""
    sum_anti_diag!(v, A)

    Sum the matrix elements over anti-diagonal directions, including
all super/sub ones
"""
function sum_anti_diag!(v::AbstractVector, A::AbstractMatrix)
    row, col = size(A)
    row == col || @error "Non-sqaure matrix"

    for (idx, i) = enumerate(-col + 1 : col - 1)
        if i < 0
            v[idx] = sum([A[j, col + 1 + i - j] for j = 1 : col + i])
        elseif i > 0
            v[idx] = sum([A[j, col + 1 + i - j] for j = 1 + i : col])
        else
            v[idx] =sum([A[j, col + 1 + i - j] for j = 1 : col])
        end
    end

    return v
end

"""
    sum_diag!(v, A)

    Sum the matrix elements over diagonal directions, including
all super/sub ones
"""
function sum_diag!(v::AbstractVector, A::AbstractMatrix)
    row, col = size(A)
    row == col || @error "Non-sqaure matrix"

    for (idx, i) = enumerate(-col + 1 : col - 1)
        if i < 0
            v[idx] = sum([A[j, j + i] for j = -i + 1 : col])
        elseif i > 0
            v[idx] = sum([A[j - i, j] for j = i + 1 : col])
        else
            v[idx] =sum([A[j, j] for j = 1 : col])
        end
    end

    return v
end
