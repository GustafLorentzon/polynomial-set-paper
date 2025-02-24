using LinearAlgebra

# Rowscaling

function rowscaling(HA, HB, C, k, alpha)

    m,_ = size(HA)

    HA2 = copy(HA);
    HB2 = copy(HB);
    C2  = copy(C);

    # Scale row
    HA2[k,:] = HA[k,:]*alpha;

    # Scale Q_k+2
    if k!=m
        HA2[k+1:end,k+2] = HA[k+1:end,k+2] / alpha;
        HB2[k+1:end,k+2] = HB[k+1:end,k+2] / alpha;
    end
    C2[k+2] = C[k+2] / alpha;

    return HA2, HB2, C2;
end

function col1shift(HA, HB, C, k::Int64, alpha)
    m,_ = size(HA)
    HA2 = copy(HA);
    HB2 = copy(HB);
    C2  = copy(C);

    # Shift identity element
    HA2[k,1] = HA[k,1] + alpha;


    # Compensate
    for j = 1:k+1
        for i = k+1:m
            HA2[i,j] = HA[i,j] - alpha * HA[i,k+2] * HB[k,j];
            HB2[i,j] = HB[i,j] - alpha * HB[i,k+2] * HB[k,j];
        end
        C2[j] = C[j] - alpha * C[k+2] * HB[k,j];
    end

    return HA2, HB2, C2;
end

function shift_ab22(HA, HB, C, alpha)
    HA2, HB2, C2 = copy(HA), copy(HB), copy(C);
    
    
    # shift a22 and b22 elements
    HA2[2,2] = HA[2,2] + alpha 
    HB2[2,2] = HB[2,2] - alpha
    
    # Compute s
    s = - alpha^2 + alpha*(HB[2,2] - HA[2,2])

    # Compute updated matrices
    HA2[3:end,3] = HA[3:end,3] - s * HA[3:end,4]
    HB2[3:end,3] = HB[3:end,3] - s * HB[3:end,4]
    C2[3]        = C2[3]       - s * C2[4];

    return HA2, HB2, C2;
end



###########################################################
# DISCLAIMER: for shift_row3 to work, alpha and beta must satisfy Z equation
###########################################################
function shift_row3(HA, HB, C, alpha, beta)
    T = eltype(HA)
    m,_ = size(HA)

    HA2, HB2, C2 = copy(HA), copy(HB), copy(C);

    Z = Zcheck(HA,HB,alpha,beta);
    Y = structure_check(HA,HB)
    if abs(Z) > 1e2*eps(real(T))
        println("Invalid transformation, alpha and beta do not satisfy conditions: ")
        @printf "%.3e" Z
    end
    if abs(Y) > 1e2*eps(real(T))
        print("Invalid transformation, HA and HB do not satisfy the structure conditions: ")
        @printf "%.3e" Y
        println()
    end

    # compute s1 s2
    s1 = alpha*(HB[3,2] - HA[3,2]) - alpha^2
    s2 = beta *(HB[3,3] - HA[3,3]) - beta^2

    # Row 3 changes
    HA2[3,2] = HA[3,2] + alpha 
    HA2[3,3] = HA[3,3] + beta 
    HB2[3,2] = HB[3,2] - alpha 
    HB2[3,3] = HB[3,3] - beta 

    # Compute changes in below rows
    if m >= 4
        HA2[4:end,3] = HA[4:end,3] - s1 * HA[4:end,5]
        HB2[4:end,3] = HB[4:end,3] - s1 * HB[4:end,5]
        HA2[4:end,4] = HA[4:end,4] - s2 * HA[4:end,5]
        HB2[4:end,4] = HB[4:end,4] - s2 * HB[4:end,5]
    end

    C2[3] = C[3] - s1 * C[5]
    C2[4] = C[4] - s2 * C[5]
    
    return HA2, HB2, C2;    

end

function Zcheck(HA, HB, alpha, beta)
    T = eltype(HA)

    # compute s1 s2
    s1 = alpha*(HB[3,2] - HA[3,2]) - alpha^2
    s2 = beta *(HB[3,3] - HA[3,3]) - beta^2

    out = alpha*(HB[3,3] - HA[3,3]) + beta*(HB[3,2] - HA[3,2]) - 2*alpha*beta - s2*HA[2,2]
    return out
end

function structure_check(HA,HB)
    col1check = maximum([maximum(HA[:,1]), maximum(HB[:,1])]);
    normcheck = maximum([ maximum(diag(HA,1) .- 1.0), maximum(diag(HB,1).-1.0) ])
    b22check  = HB[2,2] 

    Y = maximum([col1check, normcheck, b22check])
    return Y
end