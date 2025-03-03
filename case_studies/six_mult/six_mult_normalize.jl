# This will read some graphs (in Float64 / ComplexF64) and then try
# make them accurate in BigFloat. Saved in corresponding _bigfloat
#
# Note that you should not load HomotopyContinuation since
# it will have a conflicting Polynomial object
include("../common/graph_syst.jl")
include("../common/pretty_print.jl")
include("../common/degopt_tools.jl")

println("Polishing exp13")
g0=import_compgraph("six_mult_deg32_exp13_complex.cgr")

g=Compgraph(Complex{BigFloat},g0);
pw=0:32
a=big(13.0);
dervec= a.^pw;

#cref=cref[findall((!).(isinteger.(get_coeffs(g,cref) )))];

#
## Don't optimize the values that already we know
(Ha,Hb,c)=get_degopt_coeffs(g);
(Hb,Ha,c)=normalize_ab22(Hb,Ha,c)


beta=Hb[3,3];
# compute s1 s2
s2 = beta *(Hb[3,3] - Ha[3,3]) - beta^2

#s1 = alpha*(Hb[3,2] - Ha[3,2]) - alpha^2

alpha=  -(Hb[3,3]-Ha[3,3]-2*beta)\(beta*(Hb[3,2]-Ha[3,2])-s2*Ha[2,2])
(Ha,Hb,c) = shift_row3(Ha, Hb, c, alpha, beta)
#

#@show should_be_zero=abs(eval_graph(graph_degopt(Degopt(Ha,Hb,c))[1],0.3)-eval_graph(g,0.3))
#asd

(g,_)=graph_degopt(Degopt(Ha,Hb,c));



cref=get_all_cref(g)
# Exactly zero and one should be removed from optimization set
cref=cref[get_coeffs(g,cref) .!= 0]; #
cref=cref[abs.(get_coeffs(g,cref) .!= 1)]


# Remove some additional
@show size(cref)
hardcoded=[13; big(13)^32/factorial(big(32))]
for v=hardcoded
    @show v
    i=argmin(abs.(get_coeffs(g,cref).- v))
    @show cref[i]
    set_coeffs!(g,[v],[cref[i]]);
    deleteat!(cref,i);
end


@show size(cref)


sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);

eval_graph(g,0.1)



# A couple of steps of Newton
x0=get_coeffs(g,cref)

x=x0;
@show size(jac(sys,x))
for j=1:10
    global x
    J=jac(sys,x)
    F=sys(x)
    delta=J\F;
    x= x-delta
    @show norm(F)
end


(Ha,Hb,c)=get_degopt_coeffs(sys.graph);


g=sys.graph
eval_graph(g,0.1) - exp(13*0.1)
gg=Compgraph(ComplexF64,g)
eval_graph(gg,0.1) - exp(13*0.1)

export_compgraph(Compgraph(ComplexF64,g),"six_mult_deg32_exp13_complex.cgr");
