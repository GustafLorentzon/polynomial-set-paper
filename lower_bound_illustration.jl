using Random
include("common/graph_syst.jl");
include("common/pretty_print.jl");
include("rename.jl")

function setup_system(T,g0,a,poly_order,cref)
    g=Compgraph(T,g0);
    pw=0:poly_order
    dervec= a.^pw;
    dervec=convert.(T,dervec);
    scaling=1 ./dervec;
    if (all(scaling .== 1))
        scaling=convert.(T,scaling);
    end
        
    sys=GraphDervecSystem(g,dervec,scaling,cref);
    return sys
end
function cref_eq(c,creflist)
    map(x-> all(c .==  x), creflist)
end

function cref_to_index(cref::Tuple,crefselect)
    m=size(crefselect,1)
    J=zeros(Int64,m)
    for j=1:m
        #T=map( c -> (crefselect[j][1] == c[1] && crefselect[j][2] == c[2]),cref)
        J[j]=findfirst(T)
    end
    return J
end


## Setup a Degopt
Random.seed!(0);
#T=BigFloat
T=BigInt
setprecision(2000);
m=7;
Ha=ones(m,m+1); Ha=triu(Ha,1); Ha=tril(Ha,1)
Hb=ones(m,m+1); Hb=triu(Hb,1); Hb=tril(Hb,1)
Ha=Matrix{Any}(Ha); Hb=Matrix{Any}(Hb); 

c=[zeros(m+1);1]

Ha .=0;
Ha[diagind(Ha).+m] .=1 # Normalized
Ha[:,2].=1;  # Second column

Hb .=0;
view(Hb,1:m,2:(m+1)) .= 1; # Fill with ones



degopt=Degopt(convert.(T,Ha),convert.(T,Hb),convert.(T,c));
(g0,cref)=graph_degopt(degopt);
rename_to_new_notation(g0);
cref=rename_to_new_notation.(cref);

pretty_print(degopt)

#cref=cref[1:m^2];
x=get_coeffs(g0,cref);




## *****
## Compute the Jacobian 
sys=setup_system(T,g0,BigInt(1),(2^m),cref);
J=jac(sys,x);
s=svdvals(J)
@show count(s  .> 1e-17)
@show count(s  .> 1e-20)
@show count(s  .> 1e-30)
@show m^2

@show eval_graph(g0,Polynomial{T}(:x))


xx,y= get_degopt_crefs(g0)


dA=Matrix{Any}(undef,m,m+1);  dA[:].=Polynomial(0);

dB=Matrix{Any}(undef,m,m+1); dB[:].=Polynomial(0);
dc=Vector{Any}(undef,m+2); dc[:].=Polynomial(0);
DD=Diagonal(factorial.(big.(0:(2^m))));
DD=BigFloat.(DD);
for i=1:m
    for j=1:i+1
        s=(Symbol("a$i"),j)
        ii=findfirst(cref_eq(s,cref))
        tt=DD\J[:,ii];

        dA[i,j]=Polynomial(round.(BigInt,tt))

        s=(Symbol("b$i"),j)
        ii=findfirst(cref_eq(s,cref))
        tt=DD\BigFloat.(J[:,ii]);
        dB[i,j]=Polynomial(round.(BigInt,tt))
        
    end    

end
for i=1:m+2
    s=(:c,i)
    ii=findfirst(cref_eq(s,cref))
    tt=DD\J[:,ii];

    dc[i]=Polynomial(round.(BigInt,tt))
end


degvec=[vec(degree.(dA)); vec(degree.(dA-dB)); vec(degree.(dc))]

terms=Vector{Any}(undef,0);
termname=Vector{Any}(undef,0);
# Case 1
for i=1:m+2
    push!(terms,dc[i])
    push!(termname,"case 1: c$i")
end

# Case 2
for i=1:m
    for j=2:i
        push!(terms,dA[i,j]);
        push!(termname,"case 2: a$i$j")        
    end
end

#for i=1:m
#    for j=2:i
#        push!(terms,dB[i,j]);
#        push!(termname,[(Symbol("b$i"),j)])        
#    end
#end

# Case 3

for i=1:m
    for j=2:i-1
        t=dA[i,j]-dB[i,j];
        if (degree.(t)>=0)
            push!(terms,t);
            push!(termname,"case 3: a$i$j-b$i$j")
        end      
    end
end

# Case 4




for i=4:m-1
    push!(terms,dA[i,i]-dB[i,i]+2*dA[i+1,i+1]+2*(dA[i,2]-dB[i,2]))
    push!(termname,"case 4a: a$i$i-b$i$i+2*a$(i+1)$(i+1)+2*(a[$i,2]-b[$i,2])")
end

push!(terms,dA[m,m]-dB[m,m]+dc[m+1] + 2*(dA[m,2]-dB[m,2]))
push!(termname,"case 4b: a$m$m-b$m$m+c$(m+1)+2*(a$(m)2-a$(m)2)")



#
#polynomial_to_latex(dA[m,m],"a_{$m,$m}")
#polynomial_to_latex(dB[m,m],"b_{$m,$m}")
#polynomial_to_latex(dA[m,m]-dB[m,m],"x_{$m,$m}")
#polynomial_to_latex(dc[m+1],"c_{$(m+1)}")
#
#println("A - B  +C")
#polynomial_to_latex(dA[m,m]-dB[m,m]+dc[m+1],"z")
#
#
#polynomial_to_latex(dA[m,2],"a_{$m,2}")
#polynomial_to_latex(dB[m,2],"b_{$m,2}")
#polynomial_to_latex(dA[m,2]-dB[m,2],"x_{$m,2}")
#
@show length(unique(degree.(terms)))


display([termname degree.(terms) pretty_leading_terms.(terms)])
