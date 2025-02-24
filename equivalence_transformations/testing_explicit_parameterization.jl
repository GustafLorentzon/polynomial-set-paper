using LinearAlgebra
using GraphMatFun
using Polynomials
using Printf
include("explicit_parameterizations.jl")

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

HA, HB, c = P8_coeff_solver(mc);

degopt_object = Degopt(HA,HB,c)
(graph, _) = graph_degopt(degopt_object)

p1 = Polynomial(mc)
X = Polynomial([0, 1])

p2 = eval_graph(graph, X);

pdiff = p2-p1
pdiffc = coeffs(pdiff)
print("Normed difference in monomial coefficients: ")
@printf "%.3e" norm(pdiffc)