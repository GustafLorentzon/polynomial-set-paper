# Computes the degree 20 polynomial with five multiplications
#
# The output is for exp(8X) is Ha,Hb,c
# and exp(X) is Ha1,Hb1,c1

using HomotopyContinuation, GraphMatFun, LinearAlgebra
include("../common/homotopy_tools.jl")
include("../common/degopt_tools.jl")

# Input: Derivatives at zero:
#

function five_mult_deg20_degopt(target_dervec)


    # Build the evaluation scheme with the given reduction
    println("Building the system with the structure according to theory");
    m=5;
    (Ha,Hb,c,_)=build_normalized_degopt(m;reduction=[[zeros(m-2);1;1],zeros(m)])
    # This structure can always be imposed
    Ha[2,2]=0;
    (Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,3] => (Hb[3,3] +1))
    #(Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,3] => 2)
    #
    degopt=Degopt(Ha,Hb,c);
    (graph,_)=graph_degopt(degopt);

    pretty_print(Ha,Hb,c)



    # compute RHS vector
    println("Setting up system")
    @var x;
    p=eval_graph(graph,x);
    rhs=build_dervec(p,x);
    lhs=target_dervec

    @show size(rhs)
    @show size(lhs)
    expr=(rhs-lhs) ./ lhs;

    expr=expand.(expr);

    @show maximum(degree.(expr))


    @show expr[end]
    println("Using greedy elimination and imposing structure to reduce the degree")
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);

    @show degree.(expr)

    # Investigate which expr to reduce
    @show expr[argmax(degree.(expr))]  # => index: 3
    # First max degree term: B₃₋₂^2*B₅₋₅*A₃₋₂^2
    expr=expand.(subs(expr, Ha[3,2] => 0));
    (Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,2] => 0)



    @show degree.(expr)
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);
    @show degree.(expr)

    # First max degree term:  B₃₋₃^3*B₃₋₂*B₅₋₅
    expr=expand.(subs(expr, Hb[3,3] => 1));
    (Ha,Hb,c)=subs((Ha,Hb,c), Hb[3,3] => 1)
    @show degree.(expr)

    @show size(expr)
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);
    @show degree.(expr)
    @show size(expr)

    println("Solving system");
    F=System(expr);

    res=solve(F);

    v=real_solutions(res)

    T=Float64
    if (length(v) ==0)
        v=solutions(res);
        T=ComplexF64
    end

    @show norm.(v)
    @show minimum(norm.(v))

    i=argmin(norm.(v))
    x=v[i];

    @show norm(x)

    (Ha,Hb,c)=subs((Ha,Hb,c), variables(F) => x);


    @show variables(expr)[end], x[end]

    # Workaround and warning
    # for https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl/issues/615
    Ha = to_number.(Ha);
    Hb = to_number.(Hb);
    c= to_number.(c);
    if (!all((!).(typeof.(Ha) .== Expression)) ||
        !all((!).(typeof.(Hb) .== Expression)) ||
        !all((!).(typeof.(c) .== Expression))  )
        @show typeof.(Ha)
        @show typeof.(Hb)
        @show typeof.(c)
        @warn("Some entries seems to still be Expression type.")
    end


    #
    degopt0=Degopt(convert.(T,Ha),convert.(T,Hb),convert.(T,c));
    return degopt0
end
