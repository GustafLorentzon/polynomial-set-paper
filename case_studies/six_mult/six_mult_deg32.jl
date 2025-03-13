# Computes the degree d polynomial with m multiplications
#
# The output is for exp(alpha*X) is Ha,Hb,c
#
#  d=32
#  m=6
#  alpha=13
#  complex

using HomotopyContinuation, GraphMatFun, LinearAlgebra, Random
include("../common/degopt_tools.jl")
include("../common/graph_syst.jl")
include("../common/tikhonov.jl")

function init_starting_values!(g,a)

    T=eltype(g)

    Ha=[0.0       1.0              0.0              0.0             0.0           0.0    0.0
        0.0       0.3              1.0              0.0             0.0           0.0    0.0
        0.0       0.6+6.0im       -0.2+0.8im        1.0             0.0           0.0    0.0
        0.0      -4.0+4.0im      -50.0+4.0im        2.0-0.4im       1.0           0.0    0.0
        0.0      -3.0+4.0im     -100.0+200.0im     10.0-20.0im      3.0-3.0im     1.0    0.0
        0.0   -2000.0-300.0im  -8000.0+6000.0im  1000.0-2000.0im  200.0-100.0im  20.0+30.0im  1.0 ]

    Hb=[0.0     1.0     0.0      0.0    0.0   0.0   0.0
        0.0     0.0     1.0      0.0    0.0   0.0   0.0
        0.0     0.3-7.0im    0.0      1.0    0.0   0.0   0.0
        0.0     1.0+0.8im    2.0-1.0im     1.0    0.0   0.0   0.0
        0.0    20.0+7.0im  -20.0-5.0im     0.2+1.0im   1.0   0.0   0.0
        0.0   -20.0-2.0im  -70.0+100.0im  10.0-20.0im  2.0-2.0im  1.0   0.0 ]

    c=[  1.0+eps()
         a
         -300000.0 + 300000.0im
         30000.0 - 200000.0im
         6000.0 - 6000.0im
         1000.0 + 1000.0im
         200.0 - 7.0im
         convert(T,big(a)^32 / factorial(big(32)))];

    (g2,_)=graph_degopt(Degopt(Ha,Hb,c))
    cref=get_all_cref(g)
    set_coeffs!(g,get_coeffs(g2,cref),cref)
end

function setup_system(T,g0,a,poly_order)
    g=Compgraph(T,g0);
    pw=0:poly_order
    dervec= a.^pw;
    dervec=convert.(T,dervec);
    sys=GraphDervecSystem(g,dervec,1 ./dervec,cref);
    return sys
end

deg=32;
a=big(13.0); # input scaling

# Initialize
(g0,_)=graph_degopt(Degopt(1im*randn(m,m+1),1im*randn(m,m+1),1im*randn(m+2)));
init_starting_values!(g0,a)

# Determine free variables
cref=get_all_cref(g0)
II=findall((!).(get_coeffs(g0,cref) .== 0))
cref=cref[II];
II=findall((!).(get_coeffs(g0,cref) .== 1))
cref=cref[II];


# Create a  GraphDervecSystem
T=ComplexF64;
sys=setup_system(T,g0,a,deg)

x0=get_coeffs(sys.graph,cref);



@show norm(sys(x0))


lambda=0.00001;
x=x0;
d=Inf;
for j=1:20000
    global x,d, tol, lambda, rad,gx,sys;

    Fv0=sys(x)
    J=jac(sys,x);


    # Tikhonov regularized Gauss-Newton
    d=tikhonov_step(J, Fv0, lambda)

    # Armijo damping
    s=1;  count=0;
    while norm(Fv0) < norm(sys(x-s*d))    && count<10 && norm(d*s)>1
        s=s/2;
        count=count+1;
    end
    x1=x-s*d;


    if (j>500) # Only do Tikhonov in the inital stages

        # Gauss-Newton step
        d=(J\Fv0);
        # Armijo damping
        s=1;  count=0;
        while norm(Fv0) < norm(sys(x-s*d))    && count<20
            s=s/2;
            count=count+1;
        end
        x2=x-s*d;
    else
        x2=x1;
    end


    # Determine if it is better to do GN or Tikhonov GN
    ii=argmin([norm(sys(x1));norm(sys(x2))])
    if (norm(sys(x1))<norm(sys(x2)))
        x=x1;
    else
        x=x2;
    end

    condJ=Float64(cond(J))
    normx=Float64(norm(x))
    normF=Float64(norm(sys(x)))
    @show (ii,j,normF,condJ, normx)

    if (normF < 1e-8)
        @show normF
        println("DONE!!!")
        break;
    end


end

println("Polishing the solution");
# Polish the solution

T=Complex{BigFloat}
sys2=setup_system(T,sys.graph,a,deg)
x=big.(x);
@show norm(sys2(x))
for j=1:10
    x[:]=x-jac(sys2,x)\sys2(x);
    @show norm(sys2(x))
end


println("Saving solution to data/ folder")
export_compgraph(deepcopy(sys2.graph),"data/exp13_deg32.cgr");
g_compressed=Compgraph(ComplexF64,deepcopy(sys2.graph));
compress_graph!(g_compressed);
gen_code("data/exp13_deg32.m",g_compressed, lang=LangMatlab());

gen_code("data/exp13_deg32.jl",g_compressed, lang=LangJulia(),funname="exp13_deg32")
