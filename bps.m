/////////////////////////////////////////////
// Bruin - Poonen - Stoll, Prop 12.17. To  //
// get our defining polynomial from theirs //
// set Q:=-Evaluate(f,[1,y,x]).            //
/////////////////////////////////////////////

load "coleman.m";

Q:=y^3 + (-x^2 - 1)*y^2 - x^3*y + x^3 + 2*x^2 + x;
p:=3; 
N:=15;
data:=coleman_data(Q,p,N);

// finite points

P1:=set_point(0,0,data);
P2:=set_point(0,1,data); 
P3:=set_point(-3,4,data); // bad but not very bad for p = 3, 5
P4:=set_point(-1,1,data);
P5:=set_point(-1,0,data);

// infinite points

P6:=set_bad_point(0,[1,0,1],true,data);
P7:=set_bad_point(0,[1,1,1],true,data);
P8:=set_bad_point(0,[1,0,0],true,data);

// Coleman integrals over divisors P-P1

IP1P2:=coleman_integrals_on_basis(P1,P2,data:delta:=20);
IP1P3:=coleman_integrals_on_basis(P1,P3,data:delta:=20);
IP1P4:=coleman_integrals_on_basis(P1,P4,data:delta:=20);
IP1P5:=coleman_integrals_on_basis(P1,P5,data:delta:=20);
IP1P6:=coleman_integrals_on_basis(P1,P6,data:delta:=20);
IP1P7:=coleman_integrals_on_basis(P1,P7,data:delta:=20);
IP1P8:=coleman_integrals_on_basis(P1,P8,data:delta:=20);

K:=pAdicField(p,N);
M:=Matrix(3,1,Vector(K,[IP1P2[i]: i in [1..3]]));
v,_:= Kernel(M);
v1:=v.1; //xi1 is then the dot product of v1 with the 3 regular 1-forms
v2:=v.2; //xi1 is then the dot product of v2 with the 3 regular 1-forms

//computing that the Coleman integrals of xi1 vanish between P1 and P = P3,...,P8
DotProduct(v1,Vector(K,[IP1P3[i]: i in [1..3]]));
DotProduct(v1,Vector(K,[IP1P4[i]: i in [1..3]]));
DotProduct(v1,Vector(K,[IP1P5[i]: i in [1..3]]));
DotProduct(v1,Vector(K,[IP1P6[i]: i in [1..3]]));
DotProduct(v1,Vector(K,[IP1P7[i]: i in [1..3]]));
DotProduct(v1,Vector(K,[IP1P8[i]: i in [1..3]]));

//computing that the Coleman integrals of xi2 vanish between P1 and P = P3,...,P8
DotProduct(v2,Vector(K,[IP1P3[i]: i in [1..3]]));
DotProduct(v2,Vector(K,[IP1P4[i]: i in [1..3]]));
DotProduct(v2,Vector(K,[IP1P5[i]: i in [1..3]]));
DotProduct(v2,Vector(K,[IP1P6[i]: i in [1..3]]));
DotProduct(v2,Vector(K,[IP1P7[i]: i in [1..3]]));
DotProduct(v2,Vector(K,[IP1P8[i]: i in [1..3]]));
