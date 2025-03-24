using HomotopyContinuation, GraphMatFun, LinearAlgebra, Random, Plots
using LaTeXStrings
include("common/degopt_tools.jl")
include("common/graph_syst.jl")

pgfplotsx()


m=7; deg=2^m;

if m == 7
    pres1 = 800
    pres2 = 900
elseif m == 6
    pres1 = 256
    pres2 = 512
else
    pres1 = 128
    pres2 = 256
end

special_val=3.333;
T=Complex{BigFloat}
setprecision(pres1);

(Ha,Hb,c)=normalized_degopt(m;special_val=special_val)

pretty_print(Degopt(Ha,Hb,c))

(g,cref)=graph_degopt(Degopt(Ha,Hb,c));
cref=cref[findall(get_coeffs(g,cref) .== special_val)];

Random.seed!(0);
dervec=(2.0).^(0:deg);

x=rand(T, size(cref,1));

######################################################
dervec=T.(dervec);
g=Compgraph(T,g);
sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);

J1=jac(sys,x);
s1=svdvals(J1)

######################################################
setprecision(pres2)

dervec=T.(dervec);
g=Compgraph(T,g);
sys=GraphDervecSystem(g,dervec, 1 ./dervec,cref);

x2=big.(x);

J2=jac(sys,x2);
s2=svdvals(J2)

######################################################
Y1 = s1[1:end]/s1[1]
Y2 = s2[1:end]/s2[1]

p=plot(Y1; yaxis=:log, linewidth = 2, color=:red, linestyle = :solid, label = "Singular Values in Initial precision")
plot!(Y2;  yaxis=:log, linewidth = 2, color=:blue,  linestyle = :dash, dash_pattern = "on 5mm off 5mm", label = "Singular Values in Enhanced precision")
plot!([m^2;m^2], [s2[end]/s2[1];1]; yaxis=:log, color=:green, linestyle = :dash, label=L"Predicted limit for non-zero singular values: $m^2$")
plot!(legend=:bottomleft, legendfontsize=12)
title!("Singular Values of Jacobian", size=(800, 600))
plot!(tickfontsize=14)
plot!(titlefontsize=24)


plot!(p)
