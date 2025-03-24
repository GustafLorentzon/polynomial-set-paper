using LinearAlgebra


#########################################################
# Row scaling
#########################################################
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

#########################################################
# Shift first column
#########################################################

function col1shift(HA, HB, C, k::Int64, alpha)
    m,_ = size(HA)
    HA2 =  copy(HA);
    HB2 =  copy(HB);
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

#########################################################
# Shift element a_22 and b_22
#########################################################

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
    if abs(Z) > 1e2*eps(real(T))
        println("Invalid transformation, alpha and beta do not satisfy conditions: " * string(Z))
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
        HA2[4:end,4] = HA[4:end,4] - s2 * HA[4:end,5]

        HB2[4:end,3] = HB[4:end,3] - s1 * HB[4:end,5]
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



function normalize_col1(Ha,Hb,c)
    m=size(Ha,1);
    for i=1:m
        (Ha,Hb,c)=col1shift(Ha, Hb, c, i, -Ha[i,1])
        (Hb,Ha,c)=col1shift(Hb, Ha, c, i, -Hb[i,1])
    end
    return (Ha,Hb,c);
end


function normalize_superdiags(Ha,Hb,c)
    m=size(Ha,1);
    for i=1:m
        j=findlast(Ha[i,:] .!= 0);
        (Ha,Hb,c)=rowscaling(Ha, Hb, c, i, 1/Ha[i,j])
        j=findlast(Hb[i,:] .!= 0);
        (Hb,Ha,c)=rowscaling(Hb, Ha, c, i, 1/Hb[i,j])
    end
    return (Ha,Hb,c);
end


function normalize_ab22(Ha,Hb,c)
    shift_ab22(Ha, Hb, c, -Ha[2,2])
end
