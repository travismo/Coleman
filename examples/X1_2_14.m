///////////////
// X_1(2,14) //
///////////////

//> Qxy<x,y>:=PolynomialRing(RationalField(),2);
//> u:=x;
//> v:=y;
//> f:=(u^2 + u)*v^3 + (u^3 + 2*u^2 - u - 1)*v^2 + (u^3 - u^2 - 4*u - 1)*v - u^2 - u;
//> Qxy!((x^2+x)^2*Evaluate(f,[x,y/(x^2+x)]));
// -x^6 + x^5*y - 3*x^5 - 3*x^4 + x^3*y^2 - 5*x^3*y - x^3 + 2*x^2*y^2 - 5*x^2*y - x*y^2 - x*y + y^3 - y^2

load "coleman.m";
Q:=-x^6 + x^5*y - 3*x^5 - 3*x^4 + x^3*y^2 - 5*x^3*y - x^3 + 2*x^2*y^2 - 5*x^2*y - x*y^2 - x*y + y^3 - y^2;
p:=5;
N:=15;
data:=coleman_data(Q,p,N);

L,v:=effective_chabauty(data:bound:=1000,e:=30);

L;
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^11),
        b := [ 1 + O(5^15), -1 + O(5^11), -1 + O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^11),
        b := [ 1 + O(5^15), O(5^15), O(5^11) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^11),
        b := [ 1 + O(5^15), O(5^11), -1 + O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^12),
        b := [ 1 + O(5^15), O(5^15), O(5^12) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^12),
        b := [ 1 + O(5^15), 1 + O(5^12), O(5^12) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^12),
        b := [ 1 + O(5^15), O(5^12), 1 + O(5^12) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^12),
        b := [ 1 + O(5^15), -1 + O(5^12), -1 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^12),
        b := [ 1 + O(5^15), O(5^12), -1 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^12),
        b := [ 1 + O(5^15), O(5^15), O(5^12) ],
        inf := false>
]
*/

Q_points(data,1000);
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 1 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 1 + O(5^15), 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^15),
        b := [ 1 + O(5^15), 0, -1 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^15),
        b := [ 1 + O(5^15), -1 + O(5^15), -1 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^15),
        b := [ 1 + O(5^15), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), -1 + O(5^15), -1 + O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, -1 + O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0 ],
        inf := true>
]
*/

v;
/*
[
    [ 1 + O(5^13), O(5^13), O(5^13), O(5^13) ],
    [ O(5^13), 1 + O(5^13), O(5^13), O(5^13) ],
    [ O(5^13), O(5^13), 1 + O(5^13), O(5^13) ],
    [ O(5^13), O(5^13), O(5^13), 1 + O(5^13) ]
]
*/
