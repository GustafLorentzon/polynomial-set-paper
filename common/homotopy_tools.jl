using HomotopyContinuation, GraphMatFun, LinearAlgebra

function build_normalized_degopt(m;reduction=(zeros(m),zeros(m)), randvar=false)

    all_vars=Vector();
    if (!randvar)
    @var A[1:m,1:m]
    @var B[1:m,1:m]
        @var c[1:m+2]
    else
        A=randn(m,m)
        B=randn(m,m)
        c=randn(m+2);
    end

    Ha=Matrix{Number}([zeros(m) tril(ones(m,m))]);
    Hb=Matrix{Number}([zeros(m) tril(ones(m,m))]);
    cc=Vector{Number}(zeros(m+2));



    for i=2:m
        for j=2:i
            if (j+reduction[1][i] <= i)
                Ha[i,j]=A[i,j]
                push!(all_vars,A[i,j])
            else
                Ha[i,j+1]=zero(Ha[i,j+1]);
            end
            if (j+reduction[2][i] <= i)
                Hb[i,j]=B[i,j]
                push!(all_vars,B[i,j])
            else
                Hb[i,j+1]=zero(Hb[i,j+1]);
            end
        end


    end

    cc[:]=c[1:m+2];
    map(s->push!(all_vars,s),c);

    return (Ha,Hb,cc,all_vars)
end


import HomotopyContinuation.subs;
function subs(T::Tuple{X,X,Y},convert) where {X,Y}
    (Ha,Hb,y)=T
    Ha_new=map(z->(z isa Variable || z isa Expression) ? subs(z,convert) : z, Ha)
    Hb_new=map(z->(z isa Variable || z isa Expression) ? subs(z,convert) : z, Hb)
    y_new=map(z->(z isa Variable || z isa Expression) ? subs(z,convert) : z, y)

    return (Ha_new,Hb_new,y_new);
end


# Elimiante the simple relations
function greedy_eliminate(degopt,expr)
    (Ha,Hb,y)=degopt;
    modified=true;
    while modified
        modified=false;
        for j=reverse(1:size(expr,1))
            thisvars=variables(expr[j])
            if (size(thisvars,1)==1)
                res=solve(System([expr[j]]));
                real_sol=real_solutions(res);
                if (size(real_sol,1)>0)
                    x=real_sol[1]
                else
                    x=solutions(res)[1];
                end

                @show j,variables(thisvars)
                (Ha,Hb,y)=subs((Ha,Hb,y),thisvars => x);
                expr=subs(expr,thisvars => x);
                modified=true;
                deleteat!(expr,j);
            end
            if (modified)
                break
            end
        end
    end
    return (Ha,Hb,y,expr)

end

# Builds a vector of derivatives of the expression p evaluated in subsval.
function build_dervec(p,x,maxder=Inf;subsval=0)

    dervec=Vector{Any}();

    push!(dervec,subs(p,x => subsval))
    pder=p;
    dercount=0
    while (subs(pder,x=>subsval) != 0) && (dercount < maxder+1)
        pder=differentiate(pder,x);
        push!(dervec,subs(pder,x => subsval));
        dercount += 1;
    end
    dervec=dervec[1:end-1] # Remove last zero

    return dervec
end

# Returns (vars, vals) that are value pairs. degoptG is a Variable / Expression degopt and degoptF are values
function match_system(degoptG, degoptF)
    vars=[];
    vals=[];
    (Ga,Gb,Gc)=get_degopt_coeffs(degoptG);
    (Fa,Fb,Fc)=get_degopt_coeffs(degoptF);
    m=size(Ga,1);


    for j=1:3
        if (j==1)
            mat1=Ga; mat2=Fa;
        elseif (j==2)
            mat1=Gb; mat2=Fb;
        else
            mat1=Gc; mat2=Fc;
        end

        for (i,v) in enumerate(mat1)
            if (to_number(mat1[i]) != to_number(mat2[i]))
                push!(vars,mat2[i])
                push!(vals,mat1[i]);
            end
        end
    end



#    II=findall(typeof.(to_number.(vars)) .== Float64)
#    @show vars[II]
#    deleteat!(vars,II)
#    deleteat!(vals,II)
#

    return (vars,vals)

end
