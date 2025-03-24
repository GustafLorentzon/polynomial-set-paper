using LinearAlgebra
using GraphMatFun
using Polynomials
using Printf
include("../common/explicit_parameterizations.jl")

########################################################
# Testing P12
########################################################

T = BigFloat

mc = randn(T, 13)

HA, HB, c = P12_coeff_solver(mc);

degopt_object = Degopt(HA,HB,c)
(graph, _) = graph_degopt(degopt_object)

p1 = Polynomial(mc)
X = Polynomial([0, 1])

p2 = eval_graph(graph, X);

pdiff = p2-p1
pdiffc = coeffs(pdiff)
print("Normed difference in monomial coefficients: ")
@printf "%.3e" norm(pdiffc)
println()

########################################################
# Testing P8
########################################################

mc = randn(T, 9)
mc = T.([0, 0, 0, 0, 0, 0, 0, 1, 1])

HA, HB, c = P8_coeff_solver(mc);

#HB =   [0 1   0 0;
#        0 0 1 0;
#        0 17/(128) 3/8 1]

#HA =   [0 1   0 0;
#        0 0.5 1 0;
#        0 -1/(128) -5/8 1]

#c = [1, 1, 0.00103759765625, 11/64, 1]


degopt_object = Degopt(HA,HB,c)
(graph, _) = graph_degopt(degopt_object)

p1 = Polynomial(mc)
X = Polynomial([0, 1])

p2 = eval_graph(graph, X);

pdiff = p2-p1
pdiffc = coeffs(pdiff)
print("Normed difference in monomial coefficients: ")
@printf "%.3e" norm(pdiffc)