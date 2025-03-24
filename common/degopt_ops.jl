import Base.^
# Form a degopt whose output is the square of the old output.
function ^(degopt::Degopt,k::Integer)
    if (k != 2)
        error("This power is not supported");
    end

    (Ha,Hb,c)=get_degopt_coeffs(degopt);


    m=size(Ha,1);

    # Keep the same Ha,Hb  but add one row
    Ha_new=similar(Ha,size(Ha,1)+1,size(Ha,2)+1)
    Ha_new[:] .=0
    Ha_new[1:m,1:m+1]=Ha;
    Ha_new[m+1,1:(m+2)] .= c; # Row contains c


    Hb_new=similar(Hb,size(Hb,1)+1,size(Hb,2)+1)
    Hb_new[:] .=0
    Hb_new[1:m,1:m+1]=Hb;
    Hb_new[m+1,1:(m+2)] .= c; # Row contains c

    # New c is just the square of the last row
    c_new=similar(c,m+3)
    c_new[1:m+2] .=0;
    c_new[m+3]=1;

    return Degopt(Ha_new,Hb_new,c_new)

end


function scale_input(degopt::Degopt,factor)
    (Ha,Hb,c)=get_degopt_coeffs(degopt);
    Ha=copy(Ha);
    Hb=copy(Hb);
    c=copy(c);
    Ha[:,2] .*= factor;
    Hb[:,2] .*= factor;
    c[2] *= factor
    return Degopt(Ha,Hb,c)
end

function scale_output(degopt::Degopt,factor)
    (Ha,Hb,c)=get_degopt_coeffs(degopt);
    Ha=copy(Ha);
    Hb=copy(Hb);
    c=copy(c);

    c[:] .*= factor

    return Degopt(Ha,Hb,c)
end



function normalized_degopt(m;reduction=[zeros(Int64,m),zeros(Int64,m)],special_val=3.0)

    T=typeof(special_val);

    Ha=zeros(T,m,m+1);
    Hb=zeros(T,m,m+1);

    reda=Int.(reduction[1]);
    redb=Int.(reduction[2]);

    for i=1:m
        for j=2:(i-reda[i])
            Ha[i,j]=special_val
        end
        Ha[i,i-reda[i]+1]=1;
        for j=2:(i+1-redb[i])
            Hb[i,j]=special_val
        end
        Hb[i,i-redb[i]+1]=1;
    end
    c=special_val*ones(T,m+2);
    return (Ha,Hb,c)
end
