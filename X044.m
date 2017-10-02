/////////////////////////////////////////////////////
// 3. modular curve X0(44) (plane model singular!) //
/////////////////////////////////////////////////////

load "coleman.m";

Q:=y^5+12*x^2*y^3-14*x^2*y^2+(13*x^4+6*x^2)*y-(11*x^6+6*x^4+x^2);
p:=7;
N:=20;
data:=coleman_data(Q,p,N);

P1:=set_point(1,1,data); // good point
P2:=set_bad_point(0,[1,0,0,0,0],false,data); // very bad point

coleman_integrals_on_basis(P1,P2,data:delta:=100);
