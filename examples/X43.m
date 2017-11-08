/////////////////////////////////////////////////////////////////////////////////////////////////////
// X43 from http://www.birs.ca/workshops/2017/17w5065/files/explicit-moduli-banff-2017-probs-2.pdf //
/////////////////////////////////////////////////////////////////////////////////////////////////////

// R<x,y,z>:=PolynomialRing(RationalField(),3);
// Q:=3*x*y^3 + x^3*z - 6*x^2*z^2 + 3*x*z^3 + z^4;
// Q:=Evaluate(Q,[y,x,1]);
// Q;
// 3*x^3*y + y^3 - 6*y^2 + 3*y + 1

load "coleman.m";
Q:=3*x^3*y + y^3 - 6*y^2 + 3*y + 1;
p:=7;
N:=10;
data:=coleman_data(Q,p,N);
L,v:=effective_chabauty(data,1000:e:=50);

L;							// Q-points found by effective Chabauty
//[
//    rec<recformat<x, b, inf, xt, bt, index> | 
//        x := O(7^7),
//        b := [ 1 + O(7^10), O(7^10), O(7^10) ],
//        inf := true>,
//    rec<recformat<x, b, inf, xt, bt, index> | 
//        x := O(7^10),
//        b := [ 1 + O(7^10), O(7^8), -3 + O(7^10) ],
//        inf := true>
//]

Q_points(data,1000);					// Q-points found by point search
//[
//    rec<recformat<x, b, inf, xt, bt, index> | 
//        x := 0,
//        b := [ 1 + O(7^10), 0, 0 ],
//        inf := true>,
//    rec<recformat<x, b, inf, xt, bt, index> | 
//        x := 0,
//        b := [ 1 + O(7^10), 0, -3 + O(7^10) ],
//        inf := true>
//]

v;							// vanishing differentials
//[
//    [ 1 + O(7^9), O(7^9), O(7^9) ],
//    [ O(7^9), 1 + O(7^9), O(7^9) ],
//    [ O(7^9), O(7^9), 1 + O(7^9) ]
//]


