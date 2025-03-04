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

    m=5;
    (Ha,Hb,c,_)=build_normalized_degopt(m;reduction=[[zeros(m-2);1;1],zeros(m)])

    # According to theory, we can impose the following relations without loss of generality:
    Ha[2,2]=0;
    (Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,3] => (Hb[3,3] +1))

    #(Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,3] => 2)
    
    # Construct the initial graph
    degopt=Degopt(Ha,Hb,c);
    (graph,_)=graph_degopt(degopt);

    pretty_print(Ha,Hb,c)



    @var x;
    #construct polynomial
    p=eval_graph(graph,x);

    # Construct derivatives of p in X, should goal is to match target input
    rhs=build_dervec(p,x);
    lhs=target_dervec

    @show size(rhs)
    @show size(lhs)
    expr=(rhs-lhs) ./ lhs;  # Consturct multivariate polynomial system with normalization

    expr=expand.(expr);     # Simplify system

    @show maximum(degree.(expr))


    @show expr[end]

    @show size(expr)
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);    #Eliminate single variable expressions
    @show size(expr)

    @show degree.(expr)

    # Investigate which expr to reduce
    @show expr[argmax(degree.(expr))]  # => index: 3
    # First max degree term: B₃₋₂^2*B₅₋₅*A₃₋₂^2
    expr=expand.(subs(expr, Ha[3,2] => 0));                 ## // use one degree of freedom to set a_32 to zero, simplifies maxdeg expression
    (Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,2] => 0)


                                                            ## // After imposed condition we check again for single variable expressions
    @show degree.(expr)
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);
    @show degree.(expr)

    # First max degree term:  B₃₋₃^3*B₃₋₂*B₅₋₅
    expr=expand.(subs(expr, Hb[3,3] => 1));                 ## // use one degree of freedom to set b_33 = 1, this also imposes a_33 = 0 I think.
    (Ha,Hb,c)=subs((Ha,Hb,c), Hb[3,3] => 1)
    @show degree.(expr)
                                                            ## // After imposed condition we check again for single variable expressions
    @show size(expr)
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);
    @show degree.(expr)
    @show size(expr)

    # Since we have no single varable expression left, and no extra degrees of freedom, we rely on homotopy coninutation to solve the remaining system
    F=System(expr);

    res=solve(F);

    # Only care about real solutions if they exist
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

    # Substitue solution into tables
    (Ha,Hb,c)=subs((Ha,Hb,c), variables(F) => x);


    @show variables(expr)[end], x[end]

    ###################################################################################################
    # Workaround and warning
    ###################################################################################################
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
        @warn("Some entries seems to still be Expression.")
    end
    ###################################################################################################
    ###################################################################################################

    # In case we got complex solutions we convert tables to complex.
    degopt0=Degopt(convert.(T,Ha),convert.(T,Hb),convert.(T,c));
    return degopt0
end
