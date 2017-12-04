Qx<x>:=PolynomialRing(RationalField());
RankBounds(x^5 - 2*x^4 - 2*x^3 - 2*x^2 - 3*x,3);
// 1 1

load "coleman.m";
Q:=y^3 - (x^5 - 2*x^4 - 2*x^3 - 2*x^2 - 3*x);
p:=7;
N:=20;
data:=coleman_data(Q,p,N);

P1:=set_point(1,-2,data);
P2:=set_point(0,0,data);
P3:=set_point(-1,0,data);
P4:=set_point(3,0,data);
P5:=set_bad_point(0,[1,0,0],true,data);

IP1P2,N2:=coleman_integrals_on_basis(P1,P2,data:e:=50);
IP1P3,N3:=coleman_integrals_on_basis(P1,P3,data:e:=50);
IP1P4,N4:=coleman_integrals_on_basis(P1,P4,data:e:=50);
IP1P5,N5:=coleman_integrals_on_basis(P1,P5,data:e:=50);

IP1P2;
N2;
/*
(12586493*7 + O(7^10) 19221514*7 + O(7^10) -19207436*7 + O(7^10) -10636635*7 + O(7^10) 128831118 + O(7^10) 67444962 + O(7^10) -23020322 + O(7^10) 401602170*7^-1 + O(7^10))
10;
*/

// not zero, so [(P1)-(P2)] is not torsion

K:=pAdicField(p,Minimum([N2,N3,N4,N5]));
M:=Matrix(4,1,Vector(K,[IP1P2[i]: i in [1..4]]));
v,_:= Kernel(M);
v1:=v.1; // xi1 is then the dot product of v1 with the 3 regular 1-forms
v2:=v.2; // xi2 is then the dot product of v1 with the 3 regular 1-forms
v3:=v.3; // xi3 is then the dot product of v1 with the 3 regular 1-forms

v1;
// (1 + O(7^9) O(7^9) O(7^9) -18106419 + O(7^9))
v2;
// (O(7^9) 1 + O(7^9) O(7^9) 12452015 + O(7^9))
v3;
// (O(7^9) O(7^9) 1 + O(7^9) 8834289 + O(7^9))


// Computing that the Coleman integrals of xi1 vanish between P1 and P = P3,...,P5
DotProduct(v1,Vector(K,[IP1P3[i]: i in [1..4]]));
DotProduct(v1,Vector(K,[IP1P4[i]: i in [1..4]]));
DotProduct(v1,Vector(K,[IP1P5[i]: i in [1..4]]));
/*
O(7^10)
O(7^10)
O(7^10)
*/

// Computing that the Coleman integrals of xi2 vanish between P1 and P = P3,...,P5
DotProduct(v2,Vector(K,[IP1P3[i]: i in [1..4]]));
DotProduct(v2,Vector(K,[IP1P4[i]: i in [1..4]]));
DotProduct(v2,Vector(K,[IP1P5[i]: i in [1..4]]));
/*
O(7^10)
O(7^10)
O(7^10)
*/

// Computing that the Coleman integrals of xi3 vanish between P1 and P = P3,...,P5
DotProduct(v2,Vector(K,[IP1P3[i]: i in [1..4]]));
DotProduct(v2,Vector(K,[IP1P4[i]: i in [1..4]]));
DotProduct(v2,Vector(K,[IP1P5[i]: i in [1..4]]));
/*
O(7^10)
O(7^10)
O(7^10)
*/

// Carry out effective Chabauty automatically (without using the above computations)

L,v:=effective_chabauty(data:bound:=1000,e:=50);
L;
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(7^20),
        b := [ 1 + O(7^20), O(7^8), O(7^16) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(7^20),
        b := [ 1 + O(7^20), O(7^9), O(7^18) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 3 + O(7^20),
        b := [ 1 + O(7^20), O(7^9), O(7^18) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(7^20),
        b := [ 1 + O(7^20), O(7^9), O(7^18) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(7^9),
        b := [ 1 + O(7^20), -2 + O(7^9), 4 + O(7^9) ],
        inf := false>
]
*/

Q_points(data,1000);
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(7^20),
        b := [ 1 + O(7^20), -2 + O(7^20), 4 + O(7^20) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(7^20), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(7^20),
        b := [ 1 + O(7^20), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 3 + O(7^20),
        b := [ 1 + O(7^20), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(7^20), 0, 0 ],
        inf := true>
]
*/

v;
/*
[
    [ 1 + O(7^10), O(7^10), O(7^10), 22247188 + O(7^10) ],
    [ O(7^10), 1 + O(7^10), O(7^10), -27901592 + O(7^10) ],
    [ O(7^10), O(7^10), 1 + O(7^10), -71872925 + O(7^10) ]
]
*/
