using Random
include("common/graph_syst.jl");
include("common/pretty_print.jl");
include("rename.jl")

function setup_system(T,g0,a,poly_order,cref)
    g=Compgraph(T,g0);
    pw=0:poly_order
    dervec= a.^pw;
    dervec=convert.(T,dervec);
    sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);
    return sys
end
function cref_to_index(cref::Tuple,crefselect)
    m=size(crefselect,1)
    J=zeros(Int64,m)
    for j=1:m
        T=map( c -> (crefselect[j][1] == c[1] && crefselect[j][2] == c[2]),cref)
        J[j]=findfirst(T)
    end
    return J
end


## Setup a Degopt
Random.seed!(0);
#T=BigFloat
T=BigInt
setprecision(2000);
m=4;
Ha=ones(m,m+1); Ha=triu(Ha,1); Ha=tril(Ha,1)
Hb=ones(m,m+1); Hb=triu(Hb,1); Hb=tril(Hb,1)
Ha=Matrix{Any}(Ha); Hb=Matrix{Any}(Hb); 

c=ones(m+2)

Ha[1:m-1 ,1:m] = randn(m-1,m)
Hb[1:m-1 ,1:m] = randn(m-1,m)
Ha .=0;
Hb .=0;
Ha[diagind(Ha)[2:end]].= 1
#Hb[diagind(Ha)[2:end]].= -1
#Ha[2:end-1,2] .= 1;
#Hb[2:end-1,2] .= -2;


#Ha[m,(m-1):end] .=1;

#Ha[2,2]=-1;
#Ha[3,3]=1;

#Ha[1:3,1:3].=1;
#Hb[1:3,1:3].=-1;
#Ha[2,2]=0
#Hb[2,2]=-1
#Hb[2,2]=-1;
#Hb[3,3]=-1;
#Ha[3,3]=1;
#


#Ha=randn(m,m+1)
#Hb=randn(m,m+1)
Ha[:,1] .=0
Hb[:,1] .=0 
#Ha[m+1:m+1:end] .= 1 .* (1:m)
Ha[m+1:m+1:end] .= 1 
Hb[m+1:m+1:end] .= 1
#Hb[m+1:m+1:end] .= 1 .* (1:m)
Hb[m,m]=0
#Hb[m,m-1]=-1
#Ha[m,m]=0
#Ha[m-1,m-1]=1/2;

#Hb[2:end,2].=-1;

#Hb[2,2]=-1;
#Hb[3,2]=-1;
#Hb[4,2]=-1;

view(Ha,1:m,2:(m+1)) .= 0; 
view(Hb,1:m,2:(m+1)) .= 1;

Ha[diagind(Ha)[2:end]].= 0

Ha[(m+1):(m+1):end] .= 1; # Make it normalized
Hb[(m+1):(m+1):end] .= 1; # Make it normalized


Ha[:,2].=1;

#Hb[2,2]=1;
#Hb[3,2]=0;
#Ha[3,3]=0;

#Ha[m,m]=0;
#Ha[m-1,m-1]=0;
#Ha[m,2]=-1
#Hb[m,2]=0;
#Ha[m-1,2]=1
#Ha[m-2,2]=1


Ha
#Hb[end-3,2:end-4].=1;
#Ha[end-3,2:end-4].=-1;


#Hb[1:end,1].=1;
#Ha[2,2]=1;
#Ha[3,3]=1;
#Ha[4,4]=1;
#Ha[5,5]=1;
#Ha[6,6]=1;
#Hb[4,4]=-1;

#Ha[2,2]=-1;
#Ha[3,3]=1;
c[1:m+1].=0;

degopt=Degopt(convert.(T,Ha),convert.(T,Hb),convert.(T,c));
(g0,cref)=graph_degopt(degopt);
rename_to_new_notation(g0);
cref=rename_to_new_notation.(cref);

pretty_print(degopt)

#cref=cref[1:m^2];
x=get_coeffs(g0,cref);




## *****
## Compute the Jacobian 
sys=setup_system(T,g0,1,(2^m),cref);
J=jac(sys,x);
s=svdvals(J)
@show count(s  .> 1e-17)
@show count(s  .> 1e-20)
@show count(s  .> 1e-30)
@show m^2

@show eval_graph(g0,Polynomial{T}(:x))


xx,y= get_degopt_crefs(g0)



# Extract appropriate order of the cref 
cref=deepcopy(y);
for j=1:m
    map(c-> push!(cref,c), xx[j][2][2:end-1]);    
    if (j>=3)
       map(c-> push!(cref,c), xx[j][1][2:end-1]);
    end

end

for j=1:m
    
end

cref = rename_to_new_notation.(cref);


#ii=findfirst(map(t-> all((:b3,3).== t), cref))
#deleteat!(cref,ii);



## Compute Jacobian again with new ordering. 
sys=setup_system(T,g0,1,(2^m),cref);
x=get_coeffs(g0,cref)
J=jac(sys,x);


J0=deepcopy(J)
@show count(svdvals(J) .> 1e-16)

#lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
## Normalize each column 
#for j=1:size(J,2)
#    J[:,j]=J[:,j]/J[lastnz[j],j];
#end
##
#

crefd=map(c->[c], cref)

#
##new_col4_cref=[(:y, 4), (:y, 5), (:y, 6), (:Ba6, 5), (:Bb6, 5)];
new_col_cref=[];
for j=1:m-2

    push!(new_col_cref,(:c,j+3))
    #new_col4_cref=[(:c, 4), (:c, 5), (:c, 6), (:b5, 5), (:a5, 5)];
end
push!(new_col_cref,(Symbol("b$(m)"),m))
push!(new_col_cref,(Symbol("a$(m)"),m))

##new_col4_cref=[(:y, 4), (:y, 5), (:y, 6), (:y, 7), (:Ba7, 6), (:Bb7, 6)];
#new_col4_cref=[(:y, 4), (:y, 5), (:y, 6), (:y, 7), (:y,8), (:Ba8, 7), (:Bb8, 7)];
##new_col4_cref=[(:y, 4), (:y, 5), (:y, 6), (:y, 7), (:y,8), (:y,9), (:Ba9, 8), (:Bb9, 8)];

##for j=4:m+1
##    global new_col_cref
##    t=length(new_col_cref)
##    J[:,j]=J0[:,cref_to_index(cref,new_col_cref)]*[ones(t-2);-1;1]
##    crefd[j]=new_col_cref;
##    new_col_cref=new_col_cref[2:end]
##    @show new_col_cref
##end

#    #J[5:end,4] .=0;
#@show norm(J[5:end,4])
#
#new_col_cref=new_col_cref[2:end]
#t=length(new_col_cref)
#J[:,5]=J0[:,cref_to_index(cref,new_col_cref)]*[ones(t-2);-1;1]
#crefd[5]=new_col_cref;
##J[5:end,4] .=0;
#@show norm(J[6:end,5])
#
#new_col_cref=new_col_cref[2:end]
#t=length(new_col_cref)
#J[:,6]=J0[:,cref_to_index(cref,new_col_cref)]*[ones(t-2);-1;1]
#crefd[6]=new_col_cref;
##J[5:end,4] .=0;
#@show norm(J[10:end,6])
#
#
#
#new_col_cref=new_col_cref[2:end]
#t=length(new_col_cref)
#J[:,7]=J0[:,cref_to_index(cref,new_col_cref)]*[ones(t-2);-1;1]
#crefd[7]=new_col_cref;
#@show norm(J[18:end,7])
#

lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
@show lastnz

merged=[];
for i=m:-1:3
    for j=2:i
        local s1
        s1=Symbol("b$i")
        t1=(s1,j)

        s2=Symbol("a$i");
        t2=(s2,j);

        k1=findfirst(map(tt -> all(t1 .==tt[1]), crefd))
        k2=findfirst(map(tt -> all(t2 .==tt[1]), crefd))

        if (!isnothing(k1) && !isnothing(k2))
            push!(merged,k1)
            push!(merged,k2)
            
            crefd[k1]=[t1;t2];
            J[:,k1]=J[:,k1]-J[:,k2];
        end
        
       # lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
       # # Normalize each column 
       # for j=1:size(J,2)
       #     J[:,j]=J[:,j]*sign(J[lastnz[j],j]);
       # end

        if (i==j && i>3)
            @show s1,s2
            @show crefd[k1]

            if (i==m)
                t3=(:c,i+1)
            else
                t3=(Symbol("a$(j+1)"),j+1)
            end
            
                
            k3=findfirst(map(tt -> all(t3 .==tt[1]), crefd))

            @show k1,k2,k3
            ZZ=Float64.([J[:,k1] J[:,k2] J[:,k3]]);
            #@show ZZ

#            J[:,k1]=J[:,k1]-J[:,k3];
#            crefd[k1]=[t1;t2;t3];
            

        end
#
#        lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
#        # Normalize each column 
#        for j=1:size(J,2)
#            J[:,j]=J[:,j]/J[lastnz[j],j];
#        end
        
    end
end
#ii=findfirst(map(t-> all((:c,m+1) .== t), cref))
#qq=J[:,end];
#@show findlast((!).(qq.==0))
#qq=qq-J[:,ii]
#@show findlast((!).(qq.==0))
#qq=qq+J[:,end-1]
#@show findlast((!).(qq.==0))
#qq=qq+2*J[:,end-2]
#@show findlast((!).(qq.==0))
#qq=qq+2*J[:,end-3]
#@show findlast((!).(qq.==0))
#qq=qq+2*J[:,end-4];
#@show findlast((!).(qq.==0))
#
#J[:,end]= J[:,end]-J[:,ii]+J[:,end-1]+2*J[:,end-2]+2*J[:,end-3]+2*J[:,end-4]
#
#
#
#qq=J[:,end-10];
#@show findlast((!).(qq.==0));
#qq=qq-2*J[:,end-5];
#@show findlast((!).(qq.==0));
#qq=qq+J[:,end-11];
#@show findlast((!).(qq.==0));
#qq=qq+2*J[:,end-12];
#@show findlast((!).(qq.==0));
#qq=qq+2*J[:,end-13];
#@show findlast((!).(qq.==0));
#
#J[:,end-10]=qq;
#
#lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
#
#
#qq=J[:,19];
#@show findlast((!).(qq.==0));
#qq=qq-2*J[:,19+m-2];
#@show findlast((!).(qq.==0));
#qq=qq+J[:,19-1];
#@show findlast((!).(qq.==0));
#qq=qq+2*J[:,19-2];
#@show findlast((!).(qq.==0));
#
#J[:,19]=qq;
#
#lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
#
#@show lastnz
#asd





#ii=findfirst(map(t-> all((:b3,3) .== t), cref))
#jj=findfirst(map(t-> all((:a4,4) .== t), cref))
#kk=ii-1
#
#J[:,ii]=2*J[:,jj]-J[:,ii]-J[:,kk]
#crefd[ii]=[crefd[ii]...;crefd[jj]...;crefd[kk]];
#

#ii=findfirst(map(t-> all((:b5,5) .== t), cref))
#jj=findfirst(map(t-> all((:c,6) .== t), cref))
#kk=findfirst(map(t-> all((:b5,4) .== t), cref))
#pp=kk-1
#qq=pp-1;
#J[:,ii]=J[:,jj]-J[:,ii]-J[:,kk]-2*J[:,pp]-2*J[:,qq]
#
#crefd[ii]=[crefd[ii]...;crefd[jj]...;crefd[kk]...;crefd[pp]...;crefd[qq]...]
##
#
#
#ii=findfirst(map(t-> all((:b4,4) .== t), cref))
#jj=findfirst(map(t-> all((:a5,5) .== t), cref))
#kk=ii-1;
#pp=kk-1
##qq=pp-1;
#J[:,ii]=2*J[:,jj]-J[:,ii]-J[:,kk]-2*J[:,pp]
#
#crefd[ii]=[crefd[ii]...;crefd[jj]...;crefd[kk]...;crefd[pp]...]
##
##
#
#lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
#@show lastnz
#
#
#@show length(unique(lastnz))
##
#
k_last=0;
counter =0;
for s=1:100
    global J, lastnz, k, k_last, counter, crefd
    lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
    for j=1:size(J,2)
        #J[:,j]=J[:,j]/J[lastnz[j],j];
        J[:,j]=J[:,j]*sign(J[lastnz[j],j]);
    end

    @show lastnz   
    J=J[:,sortperm(lastnz)];
    crefd=crefd[sortperm(lastnz)];
    lastnz=map(j-> findlast(J[:,j] .!= 0), 1:size(J,2))
    @show lastnz

    @show length(unique(lastnz))
    @show count(svdvals(J) .> 1e-60)

    T=findall(diff(lastnz) .== 0)

    cref_lengths=map(t->length(crefd[t])+length(crefd[t+1]),T)
    @show cref_lengths
    jj=argmin(cref_lengths)
    jjj=findall(cref_lengths[jj] .== cref_lengths)

    k=T[jj[1]]
#    if (counter > 0)
#        
#    k=findfirst(diff(lastnz) .==0)
#    else
#    k=findlast(diff(lastnz) .==0)
#    end
#    asd
    if (isnothing(k) || length(unique(lastnz)) == m^2)

        println("Success!");
        @show length(lastnz)
        @show length(unique(lastnz))
        @show m^2
        break
    end
    
    @show k
    @show lastnz[k]
    @show lastnz[k+1]

#    @show findlast( (!).((J[:,k]-J[:,k+1])  .≈ 0))
#    @show findlast( (!).((J[:,k]+J[:,k+1])  .≈ 0))

    
    l=lastnz[k]

    if length(crefd[k]) < length(crefd[k+1])    
        k_keep=k;
        k_overwrite=k+1;
    else
        k_keep=k+1;
        k_overwrite=k;
    end

    @show k_keep,k_overwrite
    
    
    if (J[l,k] == J[l,k+1])
        J[:,k_overwrite]=J[:,k]-J[:,k+1]
    else        
        J[:,k_overwrite]=J[l,k+1]*J[:,k]-J[l,k]*J[:,k+1]
    end

    @show 
    @show 
    
    crefd[k_overwrite]=unique([crefd[k];crefd[k+1]]);

    @show crefd[k]



    @show J[l,k]
    #J[abs.(J) .< 1e-30] .=0

    #J[abs.(J[:,k]) .<1e-100,k] .=0;
    
    if (k == k_last)
        counter = counter +1;
    end
    if ((k == k_last && counter >4) || norm(J[:,k])==0) # remove it
        global J
        println("Removing $k")
        @show count(svdvals(J) .> 1e-60)
        @show size(J)
        J=[J[:,1:(k-1)] J[:,(k+1):end]]
        crefd=[crefd[1:k-1] crefd[k+1:end]]
        @show count(svdvals(J) .> 1e-60)
        @show size(J)        
        counter = 0;       
    end
    if (size(J,2) < m^2)
        println("Problem!")
        break
    end
    

    k_last=k;

end

#
#
#J[:,6]=J[:,6]-J[:,5]
#J[:,7]=J[:,7]-J[:,6]
#J[:,7]=J[:,7]-J[:,5]
#
#J[:,8]=J[:,9]-J[:,8]
#
