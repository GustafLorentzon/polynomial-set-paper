# (g,_)=graph_exp_native_jl([3;;]); count(values(g.operations) .== :mult)
# gen_code("/tmp/exp_native_6mult.jl",g,funname="exp_native_6mult")
# (g,_)=graph_exp_native_jl([Rational(7);;]); count(values(g.operations) .== :mult)
# gen_code("/tmp/exp_native_7mult.jl",g,funname="exp_native_7mult")

#using StableRNGs
#rng = StableRNG(0);
#randn(rng,

include("data/exp_sas_6mult_1div.jl");
include("data/exp_sas_7mult_1div.jl");
include("data/exp16_deg42_bigfloat.jl")

using LinearAlgebra, Random, BenchmarkTools

Random.seed!(0); setprecision(128); n=200
A=randn(Complex{BigFloat},n,n); A=20*A/norm(A);

alpha=16;
E1= @btime exp16_deg42_bigfloat(A/alpha);
E2= @btime exp_sas_6mult_1div(A);
E3= @btime exp_sas_7mult_1div(A);

@show norm(E1-E2)
@show norm(E1-E3)
@show norm(E2-E3)

#expAref=mapreduce(i-> (A^i)/factorial(big(i)), +, 0:100)
# Compute a reference solution with high precision
pwcoeffs=[A^0]; for j=1:60; push!(pwcoeffs,pwcoeffs[end]*A/j);end
expAref=sum(pwcoeffs)

err1=norm(expAref-E1); # Error in our method
err2=norm(expAref-E2); # Error in standard implementation version 1
err3=norm(expAref-E3); # Error in standard implementation version 2




