///////////////
// X_1(2,16) //
///////////////

load "coleman.m";
Q:=x^4*y^2 + x^3*y^3 - x^3*y + 2*x*y^3 - 2*x*y + y^4 - 2*y^2 + 1;
p:=3;
N:=15;
data:=coleman_data(Q,p,N);

L,v:=effective_chabauty(data:bound:=1000,e:=30);

L;
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(3^4),
        b := [ 1 + O(3^15), O(3^15), O(3^15), O(3^8) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(3^10),
        b := [ 1 + O(3^15), -1 + O(3^15), 1 + O(3^15), -1 + O(3^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(3^10),
        b := [ 1 + O(3^15), O(3^15), -1 + O(3^15), 1 + O(3^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(3^9),
        b := [ 1 + O(3^15), O(3^15), O(3^15), 1 + O(3^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(3^10),
        b := [ 1 + O(3^15), 1 + O(3^15), O(3^15), O(3^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(3^10),
        b := [ 1 + O(3^15), -1 + O(3^10), 2 + O(3^10), -2 + O(3^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(3^10),
        b := [ 1 + O(3^15), -1 + O(3^15), O(3^15), O(3^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(3^10),
        b := [ 1 + O(3^15), 1 + O(3^10), -2 + O(3^10), -2 + O(3^10) ],
        inf := false>
]
*/

Q_points(data,1000);
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(3^15), -1 + O(3^15), 2 + O(3^15), -2 + O(3^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(3^15), 1 + O(3^15), -2 + O(3^15), -2 + O(3^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(3^15), -1 + O(3^15), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(3^15), 1 + O(3^15), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(3^15), 0, 0, 1 + O(3^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(3^15), -1 + O(3^15), 1 + O(3^15), -1 + O(3^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(3^15), 0, -1 + O(3^15), 1 + O(3^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(3^15), 0, 0, 0 ],
        inf := true>
]
*/

v;
/*
[
    [ 1 + O(3^11), O(3^11), O(3^11), O(3^11), O(3^11) ],
    [ O(3^11), 1 + O(3^11), O(3^11), O(3^11), O(3^11) ],
    [ O(3^11), O(3^11), 1 + O(3^11), O(3^11), O(3^11) ],
    [ O(3^11), O(3^11), O(3^11), 1 + O(3^11), O(3^11) ],
    [ O(3^11), O(3^11), O(3^11), O(3^11), 1 + O(3^11) ]
]
*/
