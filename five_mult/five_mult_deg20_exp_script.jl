# Computes the degree 20 polynomial with five multiplications
#
# The output is for exp(8X) is Ha,Hb,c
# and exp(X) is Ha1,Hb1,c1

using HomotopyContinuation, GraphMatFun, LinearAlgebra
include("five_mult_deg20.jl")


pw=0:20
a=8.0;
target_dervec= a.^pw;
degopt0=five_mult_deg20_degopt(target_dervec)

(Ha,Hb,c)=get_degopt_coeffs(degopt0);

(g0,_)=graph_degopt(degopt0)


@show eval_graph(g0,0.2)-exp(a*0.2)
# Reverse the scaling
(Ha1,Hb1,c1)=get_degopt_coeffs(normalize!(scale_input(degopt0,1/a)))

(Ha1,Hb1,c1)= normalize_superdiags(Ha1,Hb1,c1) # Normalize all rows

(g1,_)=graph_degopt(Degopt(Ha1,Hb1,c1));


@show eval_graph(g1,0.2)-exp(0.2)


# export_compgraph(g1,"exp_five_mult_deg20.cgr")
# export_compgraph(g0,"exp8_five_mult_deg20.cgr")


#
#
#a=1/8;
#a=1/4
#target_dervec=Float64.(factorial.(pw).*((a).^pw))
#sum((target_dervec./factorial.(pw)) .* (x.^pw))
#degopt=five_mult_deg20_degopt(target_dervec)
#eval_graph(graph_degopt(degopt)[1],0.1)
#
