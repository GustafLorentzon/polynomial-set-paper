using LinearAlgebra
using GraphMatFun
using Random
using Printf
include("equivalence_transformations2025.jl")

################################################
# Setting up an initial evaluation scheme
################################################

# Choose Datatype
T = BigFloat

# Choose number of multiplications
m = 6

# Generate randomized tables

HA0 = rand(T, m,m+1)*1e-1 .+2
HB0 = rand(T, m,m+1)*1e-1 .+2
c0  = rand(T, m+2)*1e-1 .+2

# Make Hessenberg
for k = 1:m-1
    HA0[k,k+2:end] .= 0
    HB0[k,k+2:end] .= 0
end

# 
#for k = 1:m
#    HA0[k,k+1] = 2+rand()*1e-3
#    HB0[k,k+1] = 2+rand()*1e-3
#end

HA = copy(HA0)
HB = copy(HB0)
c  = copy(c0)

################################################
# Basic reduction
################################################

# Row scaling loop

for k = 1:m
    # column 1 shift
    global HA, HB, c = col1shift(HA, HB, c, k, -HA[k,1])
    global HA, HB, c = col1shift(HB, HA, c, k, -HB[k,1])
    # Row scaling, assumes unreduced
    global HA, HB, c = rowscaling(HA, HB, c, k, 1/HA[k,k+1])
    global HA, HB, c = rowscaling(HB, HA, c, k, 1/HB[k,k+1])
end

##################################################
# Row 2 reduction
##################################################

HA, HB, c = shift_ab22(HA, HB, c, HB[2,2])
display(HA)
display(HB)

##################################################
# Row 3 reduction
##################################################

#D = HB - HA

# Now we want a soltion on alpha so that Z(alpha,beta) == 0
# Choose beta so that, 2*beta - (b33 - a33) does not equal zero
# We set things up so that, b33 = a33 + 1, meaning we choose

# beta = (a33 + b33 + 1) / 2


##################################################
# Solving Z equation for alpha
##################################################
# Positive r gives a33 = b33 + 1
#r = T(1/2)

# Following beta always has a solution in alpha
#beta = D[3,3]/2 + r
# Determine alpha
#alpha = ( HA[2,2] * beta^2 + beta*(D[3,2] - HA[2,2]*D[3,3]) ) / (2*beta - D[3,3])

# Final transformation
#HA, HB, c = shift_row3(HA, HB, c, alpha, beta)
#HA, HB, c = shift_row3(HA, HB, c, alpha, beta+0.001)
#HA, HB, c = shift_row3(HA, HB, c, alpha, :auto)
#HA, HB, c = shift_row3(HA, HB, c, :auto, beta)
HA, HB, c = shift_row3(HA, HB, c, :auto, :auto)

##################################################
# Graph Testing
##################################################

degopt_object0 = Degopt(HA0,HB0,c0)
(graph0,crefs)=graph_degopt(degopt_object0);

degopt_object = Degopt(HA,HB,c)
(graph,crefs)=graph_degopt(degopt_object);

# Compute error

X = randn(T, 100,100)

p0 = eval_graph(graph0, X)
p = eval_graph(graph, X)

relerr = norm(p0 - p)/norm(p0);

####################################################
# DISPLAY OUTPUT
####################################################


println("______________________________________")
println()
displayTablesBool = true
if displayTablesBool
    println("A")
    display(HA)
    println("B")
    display(HB)
    println("c")
    display(transpose(c))
    println("______________________________________")
    println()
end


print("Should be "* string(2*r) * ": a33 - b33 = ")
@printf "%-3e" HA[3,3] - HB[3,3]
println()
println("______________________________________")
println()
print("relative error: ")
@printf "%.3e" relerr
println()
println("______________________________________")