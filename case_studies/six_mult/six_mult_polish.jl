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

cref=get_all_cref(g)
# Exactly zero and one should be removed from optimization set
cref=cref[get_coeffs(g,cref) .!= 0]; #
cref=cref[abs.(get_coeffs(g,cref) .!= 1)]

# Remove some additional
@show size(cref)
hardcoded=[13; big(13)^32/factorial(big(32))]
for v=hardcoded
    i=argmin(abs.(get_coeffs(g,cref) .- v))
    set_coeffs!(g,[v],[cref[i]]);
    deleteat!(cref,i);
end




sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);

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

# Make the almost real elements completely real
@show abs(eval_graph(g,big(0.001))-exp(big(0.001)*big(13)))
II=findall(abs.(imag.(get_coeffs(g,cref))) .< 1e-50)
cref1=cref[II]
set_coeffs!(g,real.(get_coeffs(g,cref1)),cref1);
eval_graph(g,0.1)-exp(0.1*13)

@show abs(eval_graph(g,big(0.001))-exp(big(0.001)*big(13)))
#export_compgraph(sys.graph,"six_mult_deg32_exp13_complex_bigfloat.cgr");
