using GraphMatFun
using Polynomials
using Distributions
using Printf
include("degopt_equivalence_transforms.jl")

a11,  a12                    = rand([-1,1])*rand(Uniform(0.9999, 1.0001), (2,1))
a21,  a22,  a23              = rand([-1,1])*rand(Uniform(1.0003, 1.0004), (3,1))
a31,  a32,  a33,  a34        = rand([-1,1])*rand(Uniform(1.0003, 1.0004), (4,1))
a41,  a42,  a43,  a44,  a45  = rand([-1,1])*rand(Uniform(1.0003, 1.0004), (5,1))

b11,  b12                   = rand([-1,1])*rand(Uniform(1.0001, 1.0002), (2,1))
b21,  b22,  b23             = rand([-1,1])*rand(Uniform(1.0003, 1.0004), (3,1))
b31,  b32,  b33,  b34       = rand([-1,1])*rand(Uniform(1.0005, 1.0006), (4,1))
b41,  b42,  b43,  b44,  b45 = rand([-1,1])*rand(Uniform(1.0007, 1.0008), (5,1))

HA = [a11  a12  0.0  0.0  0.0;
      a21  a22  a23  0.0  0.0;
      a31  a32  a33  a34  0.0;
      a41  a42  a43  a44  a45];
HB = [b11  b12  0.0  0.0  0.0;
      b21  b22  b23  0.0  0.0;
      b31  b32  b33  b34  0.0;
      b41  b42  b43  b44  b45];
y = rand(Uniform(0, 3), 6);

HA = [0  1.0  0.0  0.0  0.0;
      0  1.0  1.0  0.0  0.0;
      0  1.0  1.0  1.0  0.0;
      0  1.0  1.0  1.0  1.0];

HB = [0  1.0  0.0  0.0  0.0;
      0  0    1.0  0.0  0.0;
      0  1.0  0    1.0  0.0;
      0  1.0  1.0  1.0  1.0];
y = [1  1  1  1  1  1];

#overwrite for testing b33 script:
#HA = [0  1    0.0  0.0  0.0;
#      0  a22  1    0.0  0.0;
#      0  a32  a33  1    0.0;
#      0  a42  a43  a44  1 ];

#HB = [0  1    0.0  0.0  0.0;
#      0  0    1    0.0  0.0;
#      0  b32  b33  1    0.0;
#      0  b42  b43  b44  1 ];
#y = rand(Uniform(0, 3), 6);


(HA2,HB2,y2)=degopt_simplify(HA,HB,y)

println("mult matrices\n")
opdegsyst1 = hcat(HA, HB)
opdegsyst2 = hcat(HA2, HB2)
display(opdegsyst1)
display(opdegsyst2)
#formatted_matmul1 = [string(@sprintf("%.2e", elem)) for elem in matmul1]
#formatted_matmul2 = [string(@sprintf("%.2e", elem)) for elem in matmul2]
#display(formatted_matmul1)
#display(formatted_matmul2)


(g,_)=graph_degopt(Degopt(HA,HB,y));
(g2,_)=graph_degopt(Degopt(HA2,HB2,y2));

x1 = Polynomial([0,1])
p1 = eval_graph(g,x1)
p2 = eval_graph(g2,x1)
#display(Transpose(coeffs(p1-p2)))

eval_pt = rand([-1,1])*rand();
println(p1(eval_pt)-p2(eval_pt))



###########################################
#(gg,_)=graph_horner_degopt(rand(Uniform(0, 1), 9));
#(HHA,HHB,yy)=get_degopt_coeffs(gg)
#HHA, HHB = 2*HHA, 2*HHB;

#(HHB2,HHA2,yy2)=elemental_supdiag(HHB, HHA, yy,2)
#(HHA2,HHB2,yy2)=degopt_simplify(HHA,HHB,yy)
#(gg,_)=graph_degopt(Degopt(HHA,HHB,yy));
#(gg2,_)=graph_degopt(Degopt(HHA2,HHB2,yy2));
#println(eval_graph(gg,x1))
#println(eval_graph(gg2,x1))
#println(eval_graph(gg,8))
#println(eval_graph(gg2,8))

