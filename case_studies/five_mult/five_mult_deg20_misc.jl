# Computes the degree 20 polynomial with five multiplications
#
#
#

using HomotopyContinuation, GraphMatFun, LinearAlgebra
include("five_mult_deg20.jl")


pw=0:20
a=1
target_dervec=Float64.(factorial.(pw).*((a).^pw))
#target_dervec= a.^pw;
degopt0=five_mult_deg20_degopt(target_dervec)
(Ha,Hb,c)=get_degopt_coeffs(degopt0);
(g0,_)=graph_degopt(degopt0)


x=0.1
@show aa=sum((target_dervec./factorial.(pw)) .* (x.^pw))
@show bb=eval_graph(g0,0.1)
@show aa-bb


export_compgraph(g0,"onediv_five_mult_deg20_bigfloat.cgr")
