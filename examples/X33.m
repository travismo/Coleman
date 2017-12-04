/////////////////////////////////////////////////////////////////////////////////////////////////////
// X33 from http://www.birs.ca/workshops/2017/17w5065/files/explicit-moduli-banff-2017-probs-2.pdf //
/////////////////////////////////////////////////////////////////////////////////////////////////////

// R<x,y,z>:=PolynomialRing(RationalField(),3);
// Q:=-x^3*y+x^2*y^2-x*y^3+3*x*z^3+3*y*z^3;
// Q:=-Evaluate(Q,[y,1,x]);
// Q;
//-3*x^3*y - 3*x^3 + y^3 - y^2 + y

load "coleman.m";
Q:=-3*x^3*y - 3*x^3 + y^3 - y^2 + y;
p:=5;
N:=10;
data:=coleman_data(Q,p,N);
L,v:=effective_chabauty(data:bound:=1000,e:=30);

L; 							// Q-points found by effective Chabauty
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^6),
        b := [ 1 + O(5^10), O(5^10), O(5^10) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^10),
        b := [ 1 + O(5^10), O(5^6), 3 + O(5^10) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^7),
        b := [ 1 + O(5^10), O(5^21), O(5^42) ],
        inf := false>
]
*/

Q_points(data,1000);					// Q-points found by point search
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^10), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^10), 0, 0 ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^10), 0, 3 + O(5^10) ],
        inf := true>
]
*/

v; 							// vanishing differentials
/*
[
    [ 1 + O(5^8), O(5^8), O(5^8) ],
    [ O(5^8), 1 + O(5^8), O(5^8) ],
    [ O(5^8), O(5^8), 1 + O(5^8) ]
]
*/
