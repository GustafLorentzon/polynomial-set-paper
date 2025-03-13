using Polynomials, GraphMatFun, LinearAlgebra
# Represents the system
# f1=(graph(0)-rhs[1])*scaling
# f2=(graph'(0)-rhs[2])*scaling
# ...
struct GraphDervecSystem
    graph
    rhs
    scaling
    crefs
end


# Evaluate the system
# c=coefficient parameters
function (sys::GraphDervecSystem)(c)
    set_coeffs!(sys.graph,c,sys.crefs);
    T=eltype(sys.graph);
    x=Polynomial{T}(:x)
    p=eval_graph(sys.graph,x)
    s=size(sys.rhs,1);
    pcoeffs=p.coeffs;
    padsize=s - size(pcoeffs,1);
    pcoeffs=[pcoeffs;zeros(padsize)]
    lhs=pcoeffs  .* convert.(T,factorial.(big.(0:s-1)))
    return (lhs-sys.rhs).*sys.scaling
end

function jac(sys::GraphDervecSystem,c)
    crefs=sys.crefs
    set_coeffs!(sys.graph,c,crefs);
    T=eltype(sys.graph);
    # Compute mixed derivatives with help of eval_jac and Polynomial
    x=Polynomial{T}(:x)
    JJ_poly=eval_jac(sys.graph,[x],cref)[1,:] # Get the polynomial Jacobian
    s=size(sys.rhs,1);
    JJ=zeros(T,s,size(crefs,1));
    for j=1:size(crefs,1)
        J0=JJ_poly[j] # Contains derivative wrt to crefs[j]
        s0=length(J0.coeffs);
        JJ[1:s0,j]=J0.coeffs.*convert.(T,factorial.(big.(0:s0-1))) .*sys.scaling[1:s0];
    end
    return JJ
end


function jac_num(sys::GraphDervecSystem,c;delta=sqrt(eps()))
    n=size(sys.rhs,1);
    m=size(c,1);
    J=zeros(eltype(c),n,m);


    for k=1:m
        ek=zeros(m);
        ek[k]=1;
        J[:,k]=(sys(x+delta*ek)-sys(x-delta*ek))/(2*delta)

    end
    return J;
end

# theta=0 => sys1
function interpolate_graph_sys(sys1,sys2,theta)

    return GraphDervecSystem(sys1.graph,
                             sys2.rhs*theta+sys1.rhs*(1-theta),
                             sys1.scaling,
                             sys1.crefs)
end




#
