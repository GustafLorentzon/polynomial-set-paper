using GraphMatFun
using LinearAlgebra
using Symbolics
using Printf
include("./equiv_trans_2025/explicit_parameterizations.jl")

########################################################
# Testing P8 inverse expansion in P8
########################################################
#T = Rational{Int64}

#mc = randn(T, 9)
#mc = T.([1, 1, 1//factorial(2), 1//factorial(3), 1//factorial(4), 1//factorial(5), 1//factorial(6), 1//factorial(7), 1//factorial(8)])

# Lets try with epsilon
@variables ϵ
@variables X

T = typeof(sum(BigFloat(0.0) + ϵ))
mc = T.([0, 0, 0, 0, 0, 0, 0, 1, ϵ])
p1 = X^7 + ϵ * X^8


################## CONSTUCT TABLES USING SOLVER ######################

HA, HB, c = P8_coeff_solver(mc);
HA = simplify.(expand.(simplify.(HA)))
HB = simplify.(expand.(simplify.(HB)))
c = simplify.(expand.(simplify.(c)))

degopt_object = Degopt(HA,HB,c)
(graph, _) = graph_degopt(degopt_object)


################## CONSTUCT OUTPUT POLYNOMIAL ######################
p2 = eval_graph(graph, X);

########### test table error for random input X and epsilon:
for i=1:1
    inx, ineps = randn(2);
    println("Randomized input:")
    print(inx); print(" "); println(ineps);
    ptemp1 = substitute(p1, Dict([ϵ => ineps, X=>inx]))
    ptemp2 = substitute(p2, Dict([ϵ => ineps, X=>inx]))
    pdiff = ptemp1 - ptemp2
    println("error:")
    println(pdiff)
end
#################### PRINTS ######################################


println("Input Tables:")
display(HA)
println()
display(HB)
println()
display(c')

println()
print("Output polynomial: ")
println(simplify(p2))
print("(1/562949953421312)*(5.62949953421312e14) = ")
println((1/562949953421312)*(5.62949953421312e14))