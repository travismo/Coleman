/////////////////////////////////////////////////////////////
// Examples from https://arxiv.org/abs/math/0508174        //
// "Twists of X(7) and primitive solutions to x^2+y^3=z^7" //
/////////////////////////////////////////////////////////////


/////////
// C1: //
/////////

load "coleman.m";
Q:=6*x^3*y+y^3+x;
p:=5;
N:=15;
data:=coleman_data(Q,p,N);
Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=25);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^4),
        b := [ 1 + O(5^15), O(5^15), O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^8),
        b := [ 1 + O(5^15), O(5^4), -6 + O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^9),
        b := [ 1 + O(5^15), O(5^3), O(5^6) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -312 + O(5^4),
        b := [ 1 + O(5^15), 1562 + O(5^5), -781 + O(5^5) ],
        inf := false>
]
> 
> Qpoints;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -15258789062 + O(5^15),
        b := [ 1 + O(5^15), 15258789062 + O(5^15), -7629394531 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0 ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, -6 + O(5^15) ],
        inf := true>
]
> 
> v;
[
    [ 1 + O(5^5), O(5^5), -868 + O(5^5) ],
    [ O(5^5), 1 + O(5^5), 16 + O(5^5) ]
]
*/

/////////
// C2: //
/////////

load "coleman.m";
Q:=3*x^3*y+y^3+2*x;
p:=5;
N:=15;
data:=coleman_data(Q,p,N);
Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=25);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^4),
        b := [ 1 + O(5^15), O(5^15), O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^8),
        b := [ 1 + O(5^15), O(5^4), -3 + O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^9),
        b := [ 1 + O(5^15), O(5^3), O(5^6) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^8),
        b := [ 1 + O(5^15), -1 + O(5^4), 1 + O(5^4) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^4),
        b := [ 1 + O(5^15), 2 + O(5^5), 4 + O(5^5) ],
        inf := false>
]
> 
> Qpoints;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^15),
        b := [ 1 + O(5^15), 2 + O(5^15), 4 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(5^15),
        b := [ 1 + O(5^15), -1 + O(5^15), 1 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0 ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, -3 + O(5^15) ],
        inf := true>
]
> 
> v;
[
    [ 1 + O(5^5), O(5^5), -1303 + O(5^5) ],
    [ O(5^5), 1 + O(5^5), -42*5^2 + O(5^5) ]
]
> 
*/


////////
// C3 //
////////

load "coleman.m";
Q:=y^3 + 2*x^3*y + 3*x; // take x:=1, y:=x, z:=y. 
p:=5;
N:=15;
data:=coleman_data(Q,p,N);

Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=25);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^4),
        b := [ 1 + O(5^15), O(5^15), O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^8),
        b := [ 1 + O(5^15), O(5^4), -2 + O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^9),
        b := [ 1 + O(5^15), O(5^3), O(5^6) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(5^5),
        b := [ 1 + O(5^15), -1 + O(5^4), 1 + O(5^4) ],
        inf := false>
]
> Qpoints;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(5^15),
        b := [ 1 + O(5^15), -1 + O(5^15), 1 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0 ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, -2 + O(5^15) ],
        inf := true>
]
> v;
[
    [ 1 + O(5^5), O(5^5), 401 + O(5^5) ],
    [ O(5^5), 1 + O(5^5), -1313 + O(5^5) ]
]
*/


////////
// C4 //
////////

load "coleman.m";
Q:=-(7*x^3*y+3*x^2-3*x*y^2+y-y^4); // take x:=x, y:=1, z:=y and multiply by -1. 
p:=5;
N:=15;
data:=coleman_data(Q,p,N);

Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=25);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), O(5^15), O(5^15), O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -67073887*5 + O(5^13),
        b := [ 1 + O(5^15), 206040388 + O(5^13), -114371956 + O(5^13), -491340178 + O(5^13) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), O(5^26), O(5^52), O(5^78) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), 1 + O(5^13), 1 + O(5^13), 1 + O(5^13) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -335919187 + O(5^13),
        b := [ 1 + O(5^15), -56720952 + O(5^13), -449916821 + O(5^13), -285425783 + O(5^13) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 605718033 + O(5^13),
        b := [ 1 + O(5^15), -306498548 + O(5^13), 62827054 + O(5^13), 541069908 + O(5^13) ],
        inf := false>
]
> Qpoints;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 1 + O(5^15), 1 + O(5^15), 1 + O(5^15) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0, 0 ],
        inf := true>
]
> v;
[
    [ 1 + O(5^14), 2394632966 + O(5^14), O(5^14) ]
]
*/

// 3 points left to rule out, same or worse for other primes, MW sieving?


////////
// C6 //
////////

load "coleman.m";
Q:=y^4+2*y^3*x+3*y^2*x^2+2*y*x^3-9*y*x^2+3*y*x-y+3*x^3-y; // take x:=y, y:=x, z:=1.
p:=5;
N:=15;
data:=coleman_data(Q,p,N);

Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=25);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), O(5^13), O(5^15), O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), -1 + O(5^13), 1 + O(5^13), -1 + O(5^13) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), O(5^39), O(5^78), O(5^117) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 102352389*5 + O(5^13),
        b := [ 1 + O(5^15), 256555568 + O(5^13), 194458874 + O(5^13), 35772932 + O(5^13) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 483340367 + O(5^13),
        b := [ 1 + O(5^15), 285439741 + O(5^13), -205118544 + O(5^13), 301342896 + O(5^13) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 415363618 + O(5^13),
        b := [ 1 + O(5^15), -513577579 + O(5^13), -280139384 + O(5^13), -267494289 + O(5^13) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 572843061 + O(5^13),
        b := [ 1 + O(5^15), 501166792 + O(5^13), 327399389 + O(5^13), -172266162 + O(5^13) ],
        inf := false>
]
> Qpoints;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0, 0 ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), -1 + O(5^15), 1 + O(5^15), -1 + O(5^15) ],
        inf := true>
]
> v;
[
    [ 1 + O(5^14), -1589967487 + O(5^14), 630135344 + O(5^14) ]
]
*/

// 4 points left to rule out, same or worse for other primes, MW sieving?


