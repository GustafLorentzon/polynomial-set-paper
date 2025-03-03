# Computes the degree 20 polynomial with five multiplications
#
# The output is for exp(8X) is Ha,Hb,c
# and exp(X) is Ha1,Hb1,c1

using HomotopyContinuation, GraphMatFun, LinearAlgebra
include("../common/homotopy_tools.jl")
include("../common/degopt_tools.jl")
include("newton.jl");
# Input: Derivatives at zero:
#

@var x
pw=0:32
a=13 + x -x;
target_dervec= a.^pw;
#degopt0=five_mult_deg20_degopt(target_dervec)


#function six_mult_deg32_degopt(target_dervec)

    m=6;
    (Ha,Hb,c,_)=build_normalized_degopt(m;reduction=[[0;0;0;0;0;0],[0;0;0;1;1;1]]);

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

Ha0=[ 0.0+0.0im        1.0+0.0im              0.0+0.0im            0.0+0.0im          0.0+0.0im          0.0+0.0im      0.0+0.0im
 0.0+0.0im  -0.896154-2.3281e-146im      1.0+0.0im            0.0+0.0im          0.0+0.0im          0.0+0.0im      0.0+0.0im
 0.0+0.0im   -1.81361-0.712437im      2.1417+0.420496im       1.0+0.0im          0.0+0.0im          0.0+0.0im      0.0+0.0im
 0.0+0.0im   -4.33534+3.86382im      3.77104+3.77591im    3.29755-0.33542im      1.0+0.0im          0.0+0.0im      0.0+0.0im
 0.0+0.0im   -2.60636+3.99429im      20.6787-0.0620928im  14.5337-22.0481im   2.7372-3.08488im      1.0+0.0im      0.0+0.0im
 0.0+0.0im   -2338.69-287.361im      668.929-1789.23im    1393.47-1871.98im  151.856-116.953im  19.6302+30.2959im  1.0+0.0im]

Hb0=[
 0.0+0.0im       1.0+0.0im                 0.0+0.0im           0.0+0.0im          0.0+0.0im      0.0+0.0im  0.0+0.0im
 0.0+0.0im   1.20385-5.35865e-147im        1.0+0.0im           0.0+0.0im          0.0+0.0im      0.0+0.0im  0.0+0.0im
 0.0+0.0im   2.71414+0.161636im      -0.158304+0.420496im      1.0+0.0im          0.0+0.0im      0.0+0.0im  0.0+0.0im
 0.0+0.0im   1.13988+0.826202im        2.94971-1.26149im       1.0+0.0im          0.0+0.0im      0.0+0.0im  0.0+0.0im
 0.0+0.0im   23.4647+6.56268im         35.9522-4.29267im   1.71318+1.05518im      1.0+0.0im      0.0+0.0im  0.0+0.0im
 0.0+0.0im  -16.8622-1.50903im         26.7978-10.8821im   15.4528-18.5417im  1.60359-2.23947im  1.0+0.0im  0.0+0.0im]

c0=[                                                                              1.0 + 1.425286908033917965653869663619293829118360775124111009975439212961340875882322e-145im
                                                                             13.0 + 1.017739483352568626428997920478549447026234110564682751082283018419067948851187e-144im
 42082.89884885038019767516201183845651810163512190992115406758945747301912004582 - 217738.9058720810244073947824544840604261374887939555952873824663922382496313259im
 37870.01268487445860451966849142232469558901549118111974445652668802834545895981 - 168679.1132827853295109750065940816385368604978489595115081772241206752667924216im
 5543.173128579214606702717379307534334780941306557731699781937647690204395969391 - 5950.146264258989069957106986456348923607351016838702057490101268166865262003712im
 1411.082664126872322510802813884491617754068075676761521678992803371988841506943 + 1174.794871877074304452262258901107035770555184394085222042786412109539228107813im
 171.6982303599238481766881774605599337519399843750126274799192171361261037217577 - 7.436477849678895798435117297618912057535594310496628846142613286437151490016107im
  1.68273422049889524781321279947747916097258283851755840695714191972989894726116 - 4.496172311869743767274510035898181010286686420110601598586088542127250330261953e-147im]



#degopt_ref=Degopt(Ha0,Hb0,c0);
#gg=import_compgraph("/home/jarl/jobb_synced/degopt_transformations/six_mult_deg32_complex.cgr");
gg=import_compgraph("six_mult_deg32_exp13_complex.cgr");
gg=Compgraph(ComplexF64,gg);
degopt_ref=Degopt(gg);
(Ha_ref,Hb_ref,c_ref)=get_degopt_coeffs(degopt_ref)
Ha_ref=round.(10*Ha_ref)/10;
Hb_ref=round.(10*Hb_ref)/10;
c_ref=round.(10*c_ref)/10;
#Ha_ref[4,4] +=0.1;
degopt_ref=Degopt(Ha_ref,Hb_ref,c_ref);

(vals,vars)=match_system(degopt,degopt_ref)

II=typeof.(to_number.(vars)) .== Variable;
vals=to_number.(vals[II]);
vars=vars[II];

expr=expand.(expr);

println("Analyzing system (start)")
F=System(expr);

if (!(all(variables(F) .== vars)))
    error("Variables in different order!")
end

x=vals



for j=1:20
    global x, delta
    JJ=jacobian(F,x);
    JJ=big.(JJ);
    fact=factorize(JJ);
    delta=fact\F(x)
    x=x-delta;
    x=convert.(ComplexF64,x)
    @show norm(x), norm(F(x))
end

asd


for j=1:10
    global x, delta
    JJ=jacobian(F,x);
    JJ=big.(JJ);
    #@show svdvals(big.(JJ))
    delta=pinv(JJ,rtol=1e-9)*F(x);
    #delta=JJ\F(x);
    x=x-delta
    x=convert.(ComplexF64,x)
    @show norm(x), norm(F(x))
end
asd
println("New round");

for j=1:10
    global x, delta
    JJ=jacobian(F,x);
    @show svdvals(big.(JJ))
    delta=pinv(JJ,rtol=1e-12)*F(x);
    x=x-delta
    @show norm(x), norm(F(x))
end

asd
@show norm(F(x0))
asd
x0=vals;
newton_res=newton2(F,x0);
asd


(Ha2,Hb2,c2)=subs((Ha,Hb,c), Vector{Variable}(vars) => vals)

Ha2=to_number.(Ha2);
Hb2=to_number.(Hb2);
c2=to_number.(c2);
get_degopt_coeffs(degopt_ref)[1]


#for j=1:5
#    x=rand(size(variables(F),1));
#    Z=jacobian(F,x)
#
#    @show cond(big.(Z))
#end


#
x=rand(size(variables(F),1)) +1im*rand(size(variables(F),1));


asd

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



asd
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
