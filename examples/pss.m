/////////////////////////////////////////////////////////////
// Examples from https://arxiv.org/abs/math/0508174        //
// "Twists of X(7) and primitive solutions to x^2+y^3=z^7" //
/////////////////////////////////////////////////////////////

// Effective Chabauty (alone) works for C1,C2,C3,C8,C9,C10.

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
    [ 1 + O(5^5), O(5^5), 903 + O(5^5) ],
    [ O(5^5), 1 + O(5^5), -609 + O(5^5) ]
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
    [ 1 + O(5^5), O(5^5), -1419 + O(5^5) ],
    [ O(5^5), 1 + O(5^5), -63 + O(5^5) ]
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
        x := 443514583 + O(5^13),
        b := [ 1 + O(5^15), 222483102 + O(5^13), 206792404 + O(5^13), 570550958 + O(5^13) ],
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
    [ 1 + O(5^14), -2409871481 + O(5^14), 1995128568 + O(5^14) ]
]
*/

// 4 points left to rule out, same or worse for other primes, MW sieving?


////////
// C7 //
////////

/*
R<x,y,z>:=PolynomialRing(RationalField(),3);
f:=-3*x^4-6*x^3*z+6*x^2*y^2-6*x^2*y*z+15*x^2*z^2-4*x*y^3-6*x*y*z^2-4*x*z^3+6*y^\
2*z^2-6*y*z^3;
-3^3*Evaluate(f,[y/3,x,1]);
36*x^3*y - 18*x^2*y^2 - 162*x^2 + 18*x*y^2 + 54*x*y + 162*x + y^4 + 6*y^3 - 45*y^2 + 36*y
*/

load "coleman.m";
Q:=36*x^3*y - 18*x^2*y^2 - 162*x^2 + 18*x*y^2 + 54*x*y + 162*x + y^4 + 6*y^3 - 45*y^2 + 36*y;
p:=11;
N:=15;
data:=coleman_data(Q,p,N);

Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=50);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(11^10),
        b := [ 1 + O(11^15), O(11^15), O(11^15), O(11^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 857596757*11 + O(11^10),
        b := [ 1 + O(11^15), -1721851564 + O(11^10), 2543821670 + O(11^10), -4056664907 + O(11^10) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(11^10),
        b := [ 1 + O(11^15), O(11^10), O(11^20), O(11^30) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -890867608*11 + O(11^10),
        b := [ 1 + O(11^15), -12883474068 + O(11^10), 10290984988 + O(11^10), 8814526617 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 41457701*11 + O(11^10),
        b := [ 1 + O(11^15), 11379647582 + O(11^10), 9783899128 + O(11^10), 5822523681 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 663467820*11 + O(11^10),
        b := [ 1 + O(11^15), -3771967483 + O(11^10), -8378968429 + O(11^10), 11495433067 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -9752088465 + O(11^10),
        b := [ 1 + O(11^15), -3573456770 + O(11^10), 1478249800 + O(11^10), -3954791432 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -7863951962 + O(11^10),
        b := [ 1 + O(11^15), 19315486777 + O(11^11), -41050390684 + O(11^11), -125928372484 + O(11^11) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 2700601362 + O(11^10),
        b := [ 1 + O(11^15), -3984588515 + O(11^10), -10955850760 + O(11^10), -8005063755 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 919007743 + O(11^10),
        b := [ 1 + O(11^15), 12691651024 + O(11^10), -4901888642 + O(11^10), 2733701372 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -998087069 + O(11^10),
        b := [ 1 + O(11^15), 1062249109 + O(11^10), 9079508211 + O(11^10), 9340330882 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -3167441587 + O(11^10),
        b := [ 1 + O(11^15), -1301782747 + O(11^10), 6316512566 + O(11^10), 9438957923 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -5907354369 + O(11^10),
        b := [ 1 + O(11^15), -1483592980 + O(11^10), -9465555934 + O(11^10), 482798832 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 8336159487 + O(11^10),
        b := [ 1 + O(11^15), -6660595561 + O(11^10), 2946474668 + O(11^10), -3505099398 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 9467011700 + O(11^10),
        b := [ 1 + O(11^15), -11919575782 + O(11^10), 4180389427 + O(11^10), -2655820946 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1142312055 + O(11^10),
        b := [ 1 + O(11^15), -4332428470 + O(11^10), 1080600196 + O(11^10), -5306342285 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(11^10),
        b := [ 1 + O(11^15), O(11^10), O(11^20), O(11^30) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -11436165630 + O(11^10),
        b := [ 1 + O(11^15), 12155899832 + O(11^10), 3299895281 + O(11^10), -12890620381 + O(11^10) ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 3364916160 + O(11^10),
        b := [ 1 + O(11^15), -36103972944 + O(11^11), 124642451134 + O(11^11), 89638763526 + O(11^11) ],
        inf := false>
]
> Qpoints;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(11^15), 0, 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(11^15),
        b := [ 1 + O(11^15), 0, 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(11^15), 0, 0, 0 ],
        inf := true>
]
> v;
[
    [ 1 + O(11^11), 107912098537 + O(11^11), 60793234130 + O(11^11) ]
]
*/

// 16 points to rule out, same or worse for other primes :( 


////////
// C8 //
////////

/*
R<x,y,z>:=PolynomialRing(RationalField(),3);
f:=2*x^4-x^3*y-12*x^2*y^2+3*x^2*z^2-5*x*y^3-6*x*y^2*z+2*x*z^3-2*y^4+6*y^3*z+3*\
y^2*z^2+2*y*z^3;
-2^3*Evaluate(f,[x,y/2,1]);
-16*x^4 + 4*x^3*y + 24*x^2*y^2 - 24*x^2 + 5*x*y^3 + 12*x*y^2 - 16*x + y^4 - 6*y^3 - 6*y^2 - 8*y
*/

load "coleman.m";
Q:=-16*x^4 + 4*x^3*y + 24*x^2*y^2 - 24*x^2 + 5*x*y^3 + 12*x*y^2 - 16*x + y^4 - 6*y^3 - 6*y^2 - 8*y;
p:=5;
N:=15;
data:=coleman_data(Q,p,N);

Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=50);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), -1 + O(5^13), 1 + O(5^13), -1 + O(5^13) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), O(5^13), O(5^26), O(5^39) ],
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
        b := [ 1 + O(5^15), -1 + O(5^15), 1 + O(5^15), -1 + O(5^15) ],
        inf := true>
]
> v;
[
    [ 1 + O(5^14), O(5^14), 1178083137 + O(5^14) ],
    [ O(5^14), 1 + O(5^14), -693820769 + O(5^14) ]
]
*/


////////
// C9 //
////////

/*
R<x,y,z>:=PolynomialRing(RationalField(),3);
f:=2*x^4+4*x^3*y-4*x^3*z-3*x^2*y^2-6*x^2*y*z+6*x^2*z^2-x*y^3-6*x*y*z^2-2*y^4+2*y^3*z-3*y^2*z^2+6*y*z^3;
(-2)^3*Evaluate(f,[x,y/(-2),1]);
-16*x^4 + 16*x^3*y + 32*x^3 + 6*x^2*y^2 - 24*x^2*y - 48*x^2 - x*y^3 - 24*x*y + y^4 + 2*y^3 + 6*y^2 + 24*y
*/

load "coleman.m";
Q:=-16*x^4 + 16*x^3*y + 32*x^3 + 6*x^2*y^2 - 24*x^2*y - 48*x^2 - x*y^3 - 24*x*y + y^4 + 2*y^3 + 6*y^2 + 24*y;
p:=5;
N:=15;
data:=coleman_data(Q,p,N);

Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=50);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^12),
        b := [ 1 + O(5^15), -2 + O(5^12), 4 + O(5^12), -8 + O(5^12) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^13),
        b := [ 1 + O(5^15), O(5^26), O(5^52), O(5^78) ],
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
        b := [ 1 + O(5^15), -2 + O(5^15), 4 + O(5^15), -8 + O(5^15) ],
        inf := true>
]
> v;
[
    [ 1 + O(5^14), O(5^14), -1436750588 + O(5^14) ],
    [ O(5^14), 1 + O(5^14), -380765633*5 + O(5^14) ]
]
*/


/////////
// C10 //
/////////

/*
R<x,y,z>:=PolynomialRing(RationalField(),3);
f:=x^3*y-x^3*z+3*x^2*z^2+3*x*y^2*z+3*x*y*z^2+3*x*z^3-y^4+y^3*z+3*y^2*z^2-12*y*z^3+3*z^4;
-Evaluate(f,[x,y,1]);
-x^3*y + x^3 - 3*x^2 - 3*x*y^2 - 3*x*y - 3*x + y^4 - y^3 - 3*y^2 + 12*y - 3
*/

load "coleman.m";
Q:=-x^3*y + x^3 - 3*x^2 - 3*x*y^2 - 3*x*y - 3*x + y^4 - y^3 - 3*y^2 + 12*y - 3;
p:=5;
N:=15;
data:=coleman_data(Q,p,N);

Qpoints:=Q_points(data,10^4);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=50);

L;
Qpoints;
v;

/*
> L;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^8),
        b := [ 1 + O(5^15), O(5^8), O(5^15), O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := O(5^9),
        b := [ 1 + O(5^15), 1 + O(5^9), 1 + O(5^9), 1 + O(5^9) ],
        inf := true>
]
> Qpoints;
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 1 + O(5^15), 1 + O(5^15), 1 + O(5^15) ],
        inf := true>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(5^15), 0, 0, 0 ],
        inf := true>
]
> v;
[
    [ 1 + O(5^12), O(5^12), -121604732 + O(5^12) ],
    [ O(5^12), 1 + O(5^12), -46644808 + O(5^12) ]
]
*/
