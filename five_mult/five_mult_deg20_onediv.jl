# Computes the degree 20 polynomial with five multiplications
# for the polynomial:
#   1+x+x^2+....x^20
#
#

using HomotopyContinuation, GraphMatFun, LinearAlgebra
include("five_mult_deg20.jl")


pw=0:20
a=1 # scaling of input
# Compute the derivatives of  p(x) in the origin
target_dervec=Float64.(factorial.(pw).*((a).^pw))
degopt0=five_mult_deg20_degopt(target_dervec)
(Ha,Hb,c)=get_degopt_coeffs(degopt0);
(g0,_)=graph_degopt(degopt0)


x=0.1
println("Verify it for x=$x");
@show aa=sum((target_dervec./factorial.(pw)) .* (x.^pw))
@show bb=eval_graph(g0,0.1)
@show aa-bb

export_compgraph(g0,"data/onediv_deg20.cgr",descr="The 20 degree Taylor approximation of 1/(1-x) with five multiplications")
