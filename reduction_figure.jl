# Generates the figure with all the possible  reductions, by combinatorial loop.
using Plots;
function get_total_deg(z1,z2) # Get the degree of the polynomial with reduction specified by z1 and z2
    m=size(z1,1);
    degvec=zeros(Int64,m+2);
    degvec[1]=0;
    degvec[2]=1;
    z1=Int.(z1);
    z2=Int.(z2);
    for j=1:m
        degvec[j+2]=degvec[j+1-z1[j]]+degvec[j+1-z2[j]];
    end
    totaldeg=degvec[end];
    return totaldeg
end

# Recursive function to find all the combinations of reductions and
# for every reduction compute the degree of the resulting polynomial
# All the found degrees are stored in data
function reduction(z1,z2,data)

    boundvec=0:(size(z1,1)-1)

    if (all(z1 .<= boundvec ) && all(z2 .<= boundvec )) # Make sure it is a valid reduction
        m=size(z1,1);
        reds=sum(z1)+sum(z2)
        e=get_total_deg(z1,z2); # Degree
        push!(data,[m reds e]);

        if (mod(size(data,1),128 )==1) # Unique data. Not too often for efficiency
            unique!(data)
            print("x");
        end


        if (sum(z1) + sum(z2) <7)  # Don't do more than 7 reductions
            z_tmp=deepcopy(z1);
            for i=1:2
                for j=2:size(z1,1)
                    if (i==1)
                        z_tmp[:]=z1;
                        z_tmp[j] += 1;
                        if (z_tmp[j] <= j) # Don't reduce more than until constant
                            reduction(z_tmp,z2,data)
                        end

                    else
                        z_tmp[:]=z2;
                        z_tmp[j] += 1;
                        if (z_tmp[j] <= j) # Don't reduce more than until constant
                            reduction(z1,z_tmp,data)
                        end
                    end
                end
            end
        end
    end

end


# Collect all reductions for all multiplications
data=[];

for m=3:7
    @show m
    println("")
    reduction(zeros(m),zeros(m),data)
end

# Plot the data
ee=0.1;
x=map(z-> z[1]+ee*z[2],data);
y=map(z-> z[3], data);
scatter(x,y)

rr=0:7
for m=3:7
    plot!(m .+ rr*ee, m^2-1 .-rr)
end
