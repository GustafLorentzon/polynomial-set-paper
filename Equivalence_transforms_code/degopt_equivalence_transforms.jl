#
# Module for simplyfying degopt matrices using equivalence transforms.
#

#
# No longer our best transformations! Sept 18 2024
#

function elemental_supdiag(HA, HB, y, row=1)
    n, m        =  size(HA);
    HA2      =  Float64.(HA); #inefficient copying, but these matrices are small anyway.
    HB2      =  Float64.(HB);

    if HA[row,row+1] !=0
        alpha         = 1/HA[row,row+1]; #inverse of superdiagonal element in HA
    else
        alpha = 1
    end
    for j = 1:row+1
        HA2[row,j] = HA[row,j]*alpha #multiply full row with alpha
    end

    if row != n
        #next row until last row:
        for i = row+1:n
            HA2[i,row+2] = HA[i,row+2]/alpha #multiply row+2 column with 1/alpha to compensate
            HB2[i,row+2] = HB[i,row+2]/alpha
        end
    end
    y2 = Float64.(y);
    #similarly in output scale with 1/alpha
    y2[row+2] = y[row+2]/alpha;

    return HA2, HB2, y2;
end

function transform_subdiag(HA,HB,y)
    n,m = size(HA);
    HA2, HB2,y2 = Float64.(HA), Float64.(HB), Float64.(y);
    for row=1:n
          HA2, HB2, y2 = elemental_supdiag(HA2, HB2, y2, row)
          HB2, HA2, y2 = elemental_supdiag(HB2, HA2, y2, row) 
    end
    return HA2, HB2, y2
end

function diffMat(HA, HB, y, row=1)
    n, m        =  size(HA);
    HAdiff      =  zeros(n,m);
    HBdiff      =  zeros(n,m);

    alpha         = -HA[row,1]; #first element of row in HA
    HAdiff[row,1] =  alpha;     #We set -alpha in the corresponding position of diffMat

    β_c  = zeros(row+1);
    for i = 1:row+1
          β_c[i] = HB[row,i]*alpha
    end
    #β_c1 = HB[1,1]*alpha;
    #β_c2 = HB[1,2]*alpha;


    if row != n
          for i = row+1:n
                for j = 1:row+1 
                      HAdiff[i,j] =- β_c[j]* HA[i,row+2]
                      HBdiff[i,j] =- β_c[j]* HB[i,row+2]
                end
          end
    end
    ydiff    = zeros(m+1); 
    for j = 1:row+1
          ydiff[j] = -y[row+2]*β_c[j]; #not 100 about this, but shouldnt be out of bounds since m+1=n+2=rowmax+2
    end

    return HAdiff, HBdiff, ydiff;
end

function col1_transform(HA,HB,y,row = 1)
    #Not memory efficient at the moment.
    HAd1, HBd1, yd1 = diffMat(HA, HB, y, row);
    HA2 = HA + HAd1;
    HB2 = HB + HBd1;
    y2  = y  + yd1;

    #NOTICE REVERSE ORDER:
    HBd2, HAd2, yd2 = diffMat(HB2, HA2, y2, row);
    HA3 = HA2 + HAd2;
    HB3 = HB2 + HBd2;
    y3  = y2  + yd2;
    return HA3, HB3, y3;
end

function full_col1_transform(HA,HB,y)
    n,m = size(HA);
    HA2, HB2 = copy(HA), copy(HB);
    y2 = copy(y)
    for row=1:n
          HA2, HB2, y2 = col1_transform(HA2, HB2, y2, row)
    end
    return HA2, HB2, y2                                    
end

#assumes other transforms have been performed
function degopt_b22_transform(HA,HB,y)
    if HB[2,3] == 0
        return HA, HB, y;
    end
    n,m = size(HA);
    HA2, HB2,y2 = Float64.(HA), Float64.(HB), Float64.(y);
    
    # store original b22, a22
    b22_hat, a22_hat = HB2[2,2], HA2[2,2];
    # same row changes:
    HB2[2,2] = 0
    HA2[2,2] = a22_hat + b22_hat;
    
    #remaining changes:
    ab22 = a22_hat*b22_hat;
    for i = 2:n
        HA2[i,3] = HA2[i,3] + HA2[i,4]*ab22;
        HB2[i,3] = HB2[i,3] + HB2[i,4]*ab22;
    end
    y2[3] = y2[3] + y2[4]*ab22;

    return HA2, HB2, y2;
end

#Next equivalence transform is to transform a32, and a33 elements

function degopt_b33_transform(HA,HB,y,beta = 0)
    alpha = HB[3,3]
    n,m = size(HA); #actually makes more sense to use m,n since m is number of mults... but keep like this for now.
    HA2, HB2,y2 = Float64.(HA), Float64.(HB), Float64.(y);

    #difference matrix for more compact formulas:
    D = HB2-HA2;

    #Computing alpha based on alpha
    beta = HB[3,3]
    #beta = 3
    #numerator = -(   beta*D[3,2] + beta*HA[2,2]*D[3,3] -beta*beta*HA2[2,2]    )
    numerator = beta*HA[2,2]*D[3,3] - beta*D[3,2] - beta*beta*HA[2,2]
    denom = (D[3,3] - 2*beta)

    alpha = numerator/denom;

    HA2[3,2]+= alpha
    HB2[3,2]-= alpha

    HA2[3,3]+= beta
    HB2[3,3]-= beta

    #Finally, the rest of the rows also need to be fixed, if this is wrong we should still only see a difference in A^3

    S1 = alpha*D[3,2] - alpha*alpha
    S2 = beta *D[3,3] - beta*beta
    for k = 4:n
        HA2[k,3] = HA[k,3] - HA[k,5]*S1
        HA2[k,4] = HA[k,4] - HA[k,5]*S2
        HB2[k,3] = HB[k,3] - HB[k,5]*S1
        HB2[k,4] = HB[k,4] - HB[k,5]*S2
    end
    y2[3] = y[3] - y[5]*S1
    y2[4] = y[4] - y[5]*S2

    return HA2, HB2, y2;
end


function degopt_simplify(HA,HB,y)
    HA2, HB2,y2 = Float64.(HA), Float64.(HB), Float64.(y);
    HA2, HB2, y2 = full_col1_transform(HA2,HB2,y2);
    HA2, HB2, y2 = transform_subdiag(HA2,HB2,y2)
    HA2, HB2, y2 = degopt_b22_transform(HA2,HB2,y2);
    HA2, HB2, y2 = degopt_b33_transform(HA2,HB2,y2);
    return HA2, HB2, y2;
end

#Finding Coeffs based on monomial
function mono_to_optdeg(c_mono)
    g0, g1, g2, g3, g4, g5, g6, g7, g8 = c_mono;
    y = zeros(eltype(c_mono), 5);
    y[1:2] = [g0 g1];
    y[5]   = g8
    a22 = g7/(2*g8)
    g3hat = g3 / g8
    g4hat = g4/g8 - (g5/g8-a22*(g6/g8 - a22*a22))*a22
    g5hat = g5/g8 - (g6/g8-a22*a22)*a22
    g6hat = g6 / g8 - a22*a22
    
    a33 = g6hat
    y[4] = y[5]*g4hat
    b32 = (g3hat-a22*g4hat)/g6hat
    a32 = g5hat - (g3hat - a22*g4hat)/g6hat
    y[3]   = g2 - y[5]*a32*b32

    Ha = [0  1   0  0;
          0 a22 1   0;
          0 a32 a33 1];

    Hb = [0 1   0   0;
          0 0   1   0;
          0 b32 0  1]

    return Ha, Hb, y;
end