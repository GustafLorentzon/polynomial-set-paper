# Computes the degree 20 polynomial with five multiplications
#
# The output is for exp(8X) is Ha,Hb,c
# and exp(X) is Ha1,Hb1,c1

using HomotopyContinuation, GraphMatFun, LinearAlgebra
include("../common/homotopy_tools.jl")
include("../common/degopt_tools.jl")
# Input: Derivatives at zero:
#

@var x
pw=0:30
a=12 + x -x;
target_dervec= a.^pw;
#degopt0=five_mult_deg20_degopt(target_dervec)


#function six_mult_deg32_degopt(target_dervec)

    m=6;
    (Ha,Hb,c,_)=build_normalized_degopt(m;reduction=[[0;0;0;0;0;0],[0;0;1;0;1;1]]);

    # According to theory:
    # Ha[2,2]=0;
    #(Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,3] => (Hb[3,3] +1))
    #(Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,3] => 2)
    #
    degopt=Degopt(Ha,Hb,c);
    (graph,_)=graph_degopt(degopt);

    pretty_print(Ha,Hb,c)



    p=eval_graph(graph,x);
    rhs=build_dervec(p,x);
    lhs=target_dervec

    @show size(rhs)
    @show size(lhs)
    expr=(rhs-lhs) ./ lhs;

    expr=expand.(expr);


    @show maximum(degree.(expr))
    @show degree.(expr)


    @show expr[end]

    @show size(expr)
    @show size(variables(expr))
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);
    @show size(expr)
    @show size(variables(expr))



    @show degree.(expr)

    # Investigate which expr to reduce
    #@show expr[argmax(degree.(expr))]  # => index: 3
    # First max degree term: B₃₋₂^2*B₅₋₅*A₃₋₂^2

    #elimvar=variables(expr[8])[10]
    #expr=expand.(subs(expr, elimvar => 0));
    #(Ha,Hb,c)=subs((Ha,Hb,c), elimvar => 0)

    expr=expand.(subs(expr, Ha[3,3] => Hb[3,3]+1));
    (Ha,Hb,c)=subs((Ha,Hb,c), Ha[3,3] => Hb[3,3]+1)


#elimvar=variables(expr[5])[10]
    elimvar=Hb[6,6]
    expr=expand.(subs(expr, elimvar => 1));
    (Ha,Hb,c)=subs((Ha,Hb,c), elimvar => 1)


    @show variables(expr[1])

    @show size(expr)
    @show size(variables(expr))
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);
    @show size(expr)
    @show size(variables(expr))


    @show variables(expr[1])



    elimvar=variables(expr[7])[10]
    expr=expand.(subs(expr, elimvar => 1));
    (Ha,Hb,c)=subs((Ha,Hb,c), elimvar => 1)

    elimvar=Hb[2,2]
    expr=expand.(subs(expr, elimvar => -Ha[2,2]));
    (Ha,Hb,c)=subs((Ha,Hb,c), elimvar => -Ha[2,2])

    elimvar=variables(expr[4])[6]
    expr=expand.(subs(expr, elimvar => 1));
    (Ha,Hb,c)=subs((Ha,Hb,c), elimvar => 1)


# First max degree term:  B₃₋₃^3*B₃₋₂*B₅₋₅
    #expr=expand.(subs(expr, Hb[3,3] => 1));
    #(Ha,Hb,c)=subs((Ha,Hb,c), Hb[3,3] => 1)
    @show size(expr)
    @show size(variables(expr))
    @show degree.(expr)
    (Ha,Hb,c,expr)=greedy_eliminate((Ha,Hb,c),expr);
    @show size(expr)
    @show size(variables(expr))
    @show degree.(expr)



    println("Analyzing system")
    F=System(expr);


for j=1:5
    x=rand(size(variables(F),1));
    Z=jacobian(F,x)

    @show cond(big.(Z))
end




    println("Solving system")

    res=solve(F, start_system=:polyhedral);


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
        @warn("Some entries seems to still be Expression.")
    end



    degopt0=Degopt(convert.(T,Ha),convert.(T,Hb),convert.(T,c));
    return degopt0
#end
