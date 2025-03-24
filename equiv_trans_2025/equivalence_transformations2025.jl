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

#########################################################
# row 3 transformation, now with checks
#########################################################

# Overloading for automatic solve of Z equation

#Fully automatic
function shift_row3(HA,HB,C, alpha::Symbol, beta::Symbol)
    T = eltype(HA)
    ################################
    # Solving Z equation for alpha
    ################################
    D = HB - HA;
    r = T(1/2)

    # Following beta always has a solution in alpha
    beta = D[3,3]/2 + r
    # Determine alpha
    alpha2 = ( HA[2,2] * beta^2 + beta*(D[3,2] - HA[2,2]*D[3,3]) ) / (2*beta - D[3,3])
    beta2 = beta
    # Final transformation
    HA, HB, C = shift_row3(HA, HB, c, alpha2, beta2)

end


#Try given input beta
function shift_row3(HA,HB,C, alpha::Symbol, beta::Number)
    T = eltype(HA)
    ################################
    # Solving Z equation for alpha
    ################################
    D = HB - HA;
    if (2*beta - D[3,3]) == 0
        @warn "No solution exists for the given beta, returning original scheme \n (all other betas are admissible solutions)"
        return HA, HB, c 
    end
    r = T(1/2)

    # Determine alpha based on given beta
    alpha2 = ( HA[2,2] * beta^2 + beta*(D[3,2] - HA[2,2]*D[3,3]) ) / (2*beta - D[3,3])

    # Final transformation
    HA, HB, C = shift_row3(HA, HB, c, alpha2, beta)

end

#Try given input alpha
function shift_row3(HA,HB,C, alpha::Number, beta::Symbol)
    ################################
    # Solving Z equation for beta
    ################################
    D = HB - HA;

    if HA[2,2] != 0 # In this case there is always a solution, but it might be complex.

        # Compute alpha based on beta
        p = (D[3,2] - 2*alpha)/HA[2,2] - D[3,3]
        q = alpha*D[3,3]/HA[2,2]
        discr = p^2 / 4 - q

        if discr < 0 && !(eltype(HA) <: Complex)
            @warn "WARNING: Imaginary solution beta"
        end
        
        sqrt_disc = sqrt(discr)
        x1 = -p / 2 + sqrt_disc
        x2 = -p / 2 - sqrt_disc
        beta2 = abs(x1) < abs(x2) ? x1 : x2 #assign smallest solution

    else

        if (2*alpha - D[3,3]) == 0
            @warn "No solution exists for the given alpha, returning original scheme \n (all other alphas are admissible solutions)"
            return HA, HB, c 
        end
        beta = alpha*D[3,3] / (2*alpha - D[3,2])

    end
    # Call main function
    HA, HB, C = shift_row3(HA, HB, c, alpha, beta2)
end

# MAIN: actual transformation of row 3

function shift_row3(HA, HB, C, alpha::Number, beta::Number)
    T = eltype(HA)
    m,_ = size(HA)

    HA2, HB2, C2 = copy(HA), copy(HB), copy(C);

    Z = Zcheck(HA,HB,alpha,beta);
    Y = structure_check(HA,HB)
    if abs(Z) > 1e2*eps(real(T))
        println("Invalid transformation, alpha and beta do not satisfy conditions: ")
        #@printf "%.3e" Z
        println(Z)
        println("returning original scheme")
        return HA, HB, C
    end
    if abs(Y) > 1e2*eps(real(T))
        print("Invalid transformation, HA and HB do not satisfy the structure conditions: ")
        #@printf "%.3e" Y
        println(Y)
        println("returning original scheme")
        return HA, HB, C
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
