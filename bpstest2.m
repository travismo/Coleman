/////////////////////////////////////////////
// Bruin - Poonen - Stoll, Prop 12.17. To  //
// get our defining polynomial from theirs //
// set Q:=-Evaluate(f,[1,y,x]).            //
/////////////////////////////////////////////

load "coleman.m";

Q:=y^3 + (-x^2 - 1)*y^2 - x^3*y + x^3 + 2*x^2 + x;
// p:=3; // problem with Magma intrinsic Roots
p:=5; 
N:=10;
Qp:=pAdicField(p,N);
data:=coleman_data(Q,p,N);

// finite points

P1:=set_point(-3,4,data); // bad but not very bad for p=5
P2:=set_point(0,1,data); 
P3:=set_point(0,0,data);
P4:=set_point(-1,1,data);
P5:=set_point(-1,0,data);

// infinite points

P6:=set_bad_point(0,[1,0,1],true,data);
P7:=set_bad_point(0,[1,1,1],true,data);
P8:=set_bad_point(0,[1,0,0],true,data);

// Coleman integrals over divisors P-P1

coleman_integrals_on_basis(P1,P2,data:delta:=100);
coleman_integrals_on_basis(P1,P3,data:delta:=100);
coleman_integrals_on_basis(P1,P4,data:delta:=100);
coleman_integrals_on_basis(P1,P5,data:delta:=100);
coleman_integrals_on_basis(P1,P6,data:delta:=100);
coleman_integrals_on_basis(P1,P7,data:delta:=100);
coleman_integrals_on_basis(P1,P8,data:delta:=100);
