////////////////////////////////////////////////////////
// Example from: https://arxiv.org/pdf/1711.06345.pdf //
////////////////////////////////////////////////////////

/*
Qx<x>:=RationalFunctionField(RationalField());
Qxy<y>:=PolynomialRing(Qx);
u:=x;
v:=y;
Q:=u^5*v^2 + 2*u^4*v^3 - u^4*v^2 - u^4*v + u^3*v^4 - 4*u^3*v^2 + u^3 + u^2*v^4 - 4*u^2*v^3 + 3*u^2*v - 2*u*v^3 + 4*u*v^2 - u + v^2 - v;
Q:=(x^3+x^2)^3*Evaluate(Q,y/(x^3+x^2));;

Q;
y^4 + (2*x^4 - 4*x^2 - 2*x)*y^3 + (x^8 - 5*x^6 - 4*x^5 + 4*x^4 + 5*x^3 + x^2)*y^2 + (-x^10 - 2*x^9 + 2*x^8 + 6*x^7 + 2*x^6 - 2*x^5 - x^4)*y + x^12 + 3*x^11 + 2*x^10 - 2*x^9 - 3*x^8 - x^7;
*/

load "coleman.m";
Q:=y^4 + (2*x^4 - 4*x^2 - 2*x)*y^3 + (x^8 - 5*x^6 - 4*x^5 + 4*x^4 + 5*x^3 + x^2)*y^2 + (-x^10 - 2*x^9 + 2*x^8 + 6*x^7 + 2*x^6 - 2*x^5 - x^4)*y + x^12 + 3*x^11 + 2*x^10 - 2*x^9 - 3*x^8 - x^7;
p:=5; 
N:=12;
data:=coleman_data(Q,p,N);

Qpoints:=Q_points(data,1000);

L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=40);

L;
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^5),
        b := [ 1 + O(5^12), O(5^10), O(5^5), O(5^5) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^4),
        b := [ 1 + O(5^12), O(5^4), -1 + O(5^4), -1 + O(5^4) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^5),
        b := [ 1 + O(5^12), 1 + O(5^5), 1 + O(5^5), 3 + O(5^5) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^4),
        b := [ 1 + O(5^12), 1 + O(5^4), O(5^8), 1 + O(5^8) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^5),
        b := [ 1 + O(5^12), O(5^5), O(5^5), -1 + O(5^5) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^5),
        b := [ 1 + O(5^12), O(5^10), O(5^10), O(5^5) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^6),
        b := [ 1 + O(5^12), O(5^3), -1 + O(5^3), 1 + O(5^3) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(5^5),
        b := [ 1 + O(5^12), O(5^5), O(5^5), O(5^5) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(5^5),
        b := [ 1 + O(5^12), 2 + O(5^5), O(5^5), 1 + O(5^6) ],
        inf := false>
]
*/

Qpoints;
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(5^12),
        b := [ 1 + O(5^12), 2 + O(5^12), 0, 1 + O(5^12) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(5^12),
        b := [ 1 + O(5^12), 0, 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^12), 1 + O(5^12), 0, 1 + O(5^12) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^12), 1 + O(5^12), 1 + O(5^12), 3 + O(5^12) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^12), 0, -1 + O(5^12), -1 + O(5^12) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^12), 0, 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^12),
        b := [ 1 + O(5^12), 0, 0, -1 + O(5^12) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^12),
        b := [ 1 + O(5^12), 0, 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^12),
        b := [ 1 + O(5^12), 0, -1 + O(5^12), 1 + O(5^12) ],
        inf := false>
]
*/

v;
/*
[
    [ 1 + O(5^7), O(5^7), O(5^7), O(5^7), 21086 + O(5^7), -904 + O(5^7) ],
    [ O(5^7), 5 + O(5^7), O(5^7), O(5^7), -18508 + O(5^7), -33491 + O(5^7) ],
    [ O(5^7), O(5^7), 1 + O(5^7), O(5^7), -38466 + O(5^7), -13196 + O(5^7) ],
    [ O(5^7), O(5^7), O(5^7), 1 + O(5^7), -4008*5 + O(5^7), 34313 + O(5^7) ]
]
*/
