include("graph_syst.jl");

dervec=sqrt.(Vector(1:5.0))
monomials=dervec ./factorial.(0:size(dervec,1)-1)
(g,cref)=graph_ps_degopt(monomials);

dervec1=round.(2*dervec)/2

#cref=cref[findall(get_coeffs(g,cref) .!=0)]


sys=GraphDervecSystem(g,dervec1,1 ./dervec1,cref);

x0=get_coeffs(g,cref)
sys(x0)

x=x0;


J=jac(sys,x)

for j=1:10
    global x

    J=jac(sys,x)
    @show get_coeffs(g,cref)
    F=sys(x)
    @show F
    x= x-J\F
    @show norm(F)
end
