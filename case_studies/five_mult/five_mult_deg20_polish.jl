# This will read some graphs (in Float64 / ComplexF64) and then try
# make them accurate in BigFloat. Saved in corresponding _bigfloat
#
# Note that you should not load HomotopyContinuation since
# it will have a conflicting Polynomial object
include("../common/graph_syst.jl")
include("../common/pretty_print.jl")

println("Polishing exp8")
g0=import_compgraph("exp8_five_mult_deg20.cgr")

g=Compgraph(BigFloat,g0);
pw=0:20
a=big(8.0);
dervec= a.^pw;

cref=get_all_cref(g)
cref=cref[findall((!).(isinteger.(get_coeffs(g,cref) )))];

# Don't optimize the values that already we know
(Ha,Hb,c)=get_degopt_coeffs(g);
Ha[4,3]=-41//16
c[end]=big(a)^20/factorial(big(20));
set_coeffs!(g,get_coeffs(graph_degopt(Degopt(Ha,Hb,c))[1]))
deleteat!(cref,findfirst(get_coeffs(g,cref) .≈ c[end]))
deleteat!(cref,findfirst(get_coeffs(g,cref) .≈ Ha[4,3]));


# Set up the system
sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);



# A couple of steps of Newton
x0=get_coeffs(g,cref)

x=x0;

for j=1:10
    global x
    J=jac(sys,x)
    F=sys(x)
    x= x-J\F
    @show norm(F)
end

degopt=Degopt(g)

g=sys.graph;

(Ha,Hb,c)=get_degopt_coeffs(g);


degopt_float64=Degopt(Compgraph(Float64,g));
latex_print_vals(degopt_float64)

export_compgraph(g,"exp8_five_mult_deg20_bigfloat.cgr")




println("Polishing onediv")
g0=import_compgraph("onediv_five_mult_deg20.cgr")

g=Compgraph(BigFloat,g0);
pw=0:20
a=big(1.0);
dervec=(factorial.(pw).*((a).^pw))

cref=get_all_cref(g)
cref=cref[findall((!).(isinteger.(get_coeffs(g,cref) )))];

# Don't optimize the values that already we know
(Ha,Hb,c)=get_degopt_coeffs(g);
Hb[2,2]=1//5;
Ha[4,3]=-27//5
set_coeffs!(g,get_coeffs(graph_degopt(Degopt(Ha,Hb,c))[1]))
deleteat!(cref,findfirst(get_coeffs(g,cref) .≈ Hb[2,2]))
deleteat!(cref,findfirst(get_coeffs(g,cref) .≈ Ha[4,3]))


# Set up the system
sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);



# A couple of steps of Newton
x0=get_coeffs(g,cref)
x=x0;
for j=1:10
    global x
    J=jac(sys,x)
    F=sys(x)
    x= x-J\F
    @show norm(F)
end

degopt=Degopt(g)

g=sys.graph;

(Ha,Hb,c)=get_degopt_coeffs(g);


degopt_float64=Degopt(Compgraph(Float64,g));
latex_print_vals(degopt_float64)

export_compgraph(g,"onediv_five_mult_deg20_bigfloat.cgr");
