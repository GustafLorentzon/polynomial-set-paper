using HomotopyContinuation, GraphMatFun, LinearAlgebra, Random
include("case_studies/common/degopt_tools.jl")
include("case_studies/common/graph_syst.jl")


m=6; deg=2^m;

special_val=3.333;
setprecision(512);

(Ha,Hb,c)=normalized_degopt(m;special_val=special_val)

pretty_print(Degopt(Ha,Hb,c))

(g,cref)=graph_degopt(Degopt(Ha,Hb,c));
cref=cref[findall(get_coeffs(g,cref) .== special_val)];

Random.seed!(0);

dervec=(2.0).^(0:deg);
x=rand(size(cref,1));

sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);
J1=jac(sys,x);


T=Complex{BigFloat}
dervec=big.(dervec);
g=Compgraph(T,g);
sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);
xb=big.(x);
Jb=jac(sys,xb);


s1=svdvals(J1)
sb=svdvals(Jb)

p=plot((s1[1:end]/s1[1]); yaxis=:log)
plot!((sb[1:end]/sb[1]); yaxis=:log)
plot!([m^2;m^2], [1e-50;1]; yaxis=:log)

plot!(p)
