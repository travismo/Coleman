/////////////////////////////////////////////
// Bruin - Poonen - Stoll, Prop 12.16. To  //
// get our defining polynomial from theirs //
// set Q:=-Evaluate(f,[x,y,1]).            //
/////////////////////////////////////////////

load "coleman.m";

Q:=y^3 + (-x^2 - x)*y^2 + x^3*y - x^2 + x;
p:=5;
N:=10;
data:=coleman_data(Q,p,N);

// finite points

P1:=set_point(1,1,data);
P2:=set_point(0,0,data);
P3:=set_point(1,0,data);

// infinite points

P4:=set_bad_point(0,[1,0,-1],true,data);
P5:=set_bad_point(0,[1,1,0],true,data);
P6:=set_bad_point(0,[1,0,0],true,data);

// Coleman integrals over divisors P-P1 (Jacobian is of rank 0, so integrals of regular 1-forms should vanish)

coleman_integrals_on_basis(P1,P2,data:delta:=100);
coleman_integrals_on_basis(P1,P3,data:delta:=100);
coleman_integrals_on_basis(P1,P4,data:delta:=100);
coleman_integrals_on_basis(P1,P5,data:delta:=100);
coleman_integrals_on_basis(P1,P6,data:delta:=100);

// try to reconstruct some of the points: 

xt,bt:=local_coord(P2,2*N,data);
x:=Evaluate(xt,p);
b:=Eltseq(Evaluate(Vector(bt),p));
P:=set_bad_point(x,b,false,data); // finite bad in residue disk of P2
torsion_points_in_disk(P1,P,data:delta:=100);

xt,bt:=local_coord(P3,2*N,data);
x:=Evaluate(xt,p);
b:=Eltseq(Evaluate(Vector(bt),p));
P:=set_bad_point(x,b,false,data); // finite bad in residue disk of P3
torsion_points_in_disk(P1,P,data:delta:=100); 


