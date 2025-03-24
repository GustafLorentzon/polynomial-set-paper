include("equivalence_transformations.jl")

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

# To impose b_33 = a_33 + 1, use: shift_row3, with two symbol inputs.