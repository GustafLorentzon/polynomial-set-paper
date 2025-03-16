# Computes the degree d polynomial with m multiplications
#
# The output is for exp(alpha*X) is Ha,Hb,c
#
#  d=42
#  m=7
#  alpha=16
#  complex

using HomotopyContinuation, GraphMatFun, LinearAlgebra, Random
include("../common/degopt_tools.jl")
include("../common/graph_syst.jl")
include("../common/tikhonov.jl")

function init_starting_values!(g,a)

    T=eltype(g)

    Ha=[  0.0    1.0              0.0           0.0                    0.0            0.0              0.0          0.0
          0.0   1.58-1.4im        1.0           0.0                    0.0            0.0              0.0          0.0
          0.0   6.38+0.891im      0.0           1.0                    0.0            0.0              0.0          0.0
          0.0  -6.81+12.5im      1.58-1.4im     4.72+0.389im           1.0            0.0              0.0          0.0
          0.0   8.25+65.1im       0.0          39.9-28.1im           4.38-1.11im      1.0              0.0          0.0
          0.0  305.0+436.0im    463.0+1930.0im     -773.0+553.0im    -87.9+96.2im    -19.4+0.0863im    1.0          0.0
          -0.0  -88.1-594.0im  46300.0-92500.0im  -13000.0+5570.0im  1090.0-4880.0im  -99.1-260.0im   -47.9-37.3im  1.0
          ];
    Hb=[ 0.0    1.0            0.0             0.0           0.0          0.0         0.0  0.0
         0.0  0.836+1.51im     1.0             0.0           0.0          0.0         0.0  0.0
         0.0  -3.35-0.188im    1.0             0.0           0.0          0.0         0.0  0.0
         0.0   6.32+0.158im   2.62-2.27im      1.0           0.0          0.0         0.0  0.0
         0.0   7.43+2.54im    19.1-1.8im      2.54+0.607im   1.0          0.0         0.0  0.0
         0.0   10.5+6.24im   -52.8-22.3im     38.2-39.2im   2.76-1.44im   1.0         0.0  0.0
         0.0  235.0+94.6im   418.0+1440.0im  156.0-242.0im  35.0+63.7im  4.03+1.08im  1.0  0.0]
    c=[1 + eps()  # + eps() in order to avoid exact singularity
       a
       2.25e+06 - 1.7e+06im
       2.3e5 + 72400.0im
       135000.0 - 60100.0im
       8760.0 + 2630.0im
       -677.0 - 1500.0im
       17.1 + 13.3im
       convert(T,big(a)^42 / factorial(big(42)))]

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

m=7;
deg=42;
a=big(16.0); # input scaling

# Initialize
(g0,_)=graph_degopt(Degopt(1im*randn(m,m+1),1im*randn(m,m+1),1im*randn(m+2)));
init_starting_values!(g0,a)
pretty_print(g0);

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


lambda=0.0001;
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
#export_compgraph(deepcopy(sys2.graph),"data/exp16_deg42.cgr");
#g_compressed=Compgraph(ComplexF64,deepcopy(sys2.graph));
#compress_graph!(g_compressed);
#gen_code("data/exp16_deg42.m",g_compressed, lang=LangMatlab());
#
#gen_code("data/exp16_deg42.jl",g_compressed, lang=LangJulia(),funname="exp16_deg42")
