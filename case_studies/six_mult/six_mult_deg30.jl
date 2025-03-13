# Computes the degree d polynomial with m multiplications
#
# The output is for exp(alpha*X) is Ha,Hb,c
#
#  d=30
#  m=6
#  alpha=30

using HomotopyContinuation, GraphMatFun, LinearAlgebra, Random
include("../common/degopt_tools.jl")
include("../common/graph_syst.jl")
include("../common/tikhonov.jl")

function init_starting_values!(g,a)

    T=eltype(g)
    Ha=[ 0.0    1.0    0.0    0.0   0.0   0.0  0.0
         0.0   -9.0    1.0    0.0   0.0   0.0  0.0
         0.0   30.0  100.0    1.0   0.0   0.0  0.0
         0.0    3.0  -10.0   -0.8   1.0   0.0  0.0
         0.0   20.0  -80.0   -0.5  -0.4   1.0  0.0
         0.0  -30.0  -30.0  -40.0  60.0  -3.0  1.0
         ];
    Hb=[ 0.0   1.0     0.0    0.0   0.0  0.0  0.0
         0.0  12.0     1.0    0.0   0.0  0.0  0.0
         0.0  -2.5     1.0    0.0   0.0  0.0  0.0
         0.0   6.0  -160.0   -2.3   1.0  0.0  0.0
         0.0  -1.0  -100.0   -1.7   1.0  0.0  0.0
         0.0  72.0  -870.0  -24.0  21.0  1.0  0.0];

    c=[1.0+1e-10
       13.0+1e-10
       240.0
       -35000.0
       50000.0
       38.0
       1300.0
       convert(T,big(a)^30 / factorial(big(30)))];



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

deg=30;
a=big(13.0); # input scaling

# Initialize
(g0,_)=graph_degopt(Degopt(randn(m,m+1),randn(m,m+1),randn(m+2)));
init_starting_values!(g0,a)

# Determine free variables
cref=get_all_cref(g0)
II=findall((!).(get_coeffs(g0,cref) .== 0))
cref=cref[II];
II=findall((!).(get_coeffs(g0,cref) .== 1))
cref=cref[II];


# Create a  GraphDervecSystem
T=Float64;
sys=setup_system(T,g0,a,deg)

x0=get_coeffs(sys.graph,cref);



@show norm(sys(x0))


lambda=0.02;
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
    while norm(Fv0) < norm(sys(x-s*d))    && count<10
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

# Polish the solution

T=BigFloat
sys2=setup_system(T,sys.graph,a,deg)
x=big.(x);
@show norm(sys2(x))
for j=1:10
    x[:]=x-jac(sys2,x)\sys2(x);
    @show norm(sys2(x))
end
