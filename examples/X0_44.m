////////////
// X0(44) //
////////////

load "coleman.m";

Q:=y^5+12*x^2*y^3-14*x^2*y^2+(13*x^4+6*x^2)*y-(11*x^6+6*x^4+x^2); // equation found by Yang
p:=7;
N:=20;
data:=coleman_data(Q,p,N);

P1:=set_point(1,1,data);
P2:=set_bad_point(0,[1,0,0,0,0],false,data);

coleman_integrals_on_basis(P1,P2,data:e:=75);
/*
(O(7^9) O(7^9) O(7^9) O(7^9) 11452342 + O(7^9) 2808875*7^-1 + O(7^9) 1314185 + O(7^9) 115220003*7^-1 + O(7^9))
9
*/

// In fact, all Coleman integrals of regular differentials between rational points seem to vanish:

Qpoints:=Q_points(data,1000);

Qpoints;
/*
[
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 1 + O(7^20),
        b := [ 1 + O(7^20), 1 + O(7^20), 1 + O(7^20), 1 + O(7^20), -2 + O(7^20) 
            ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(7^20), 0, 0, 0, 0 ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 7253842390692000 + O(7^20),
        b := [ 1 + O(7^20), -21761527172076000 + O(7^20), 15826565216055273 + 
            O(7^20), 32312570649446182 + O(7^20), -14507684781384000 + O(7^20) 
        ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -1 + O(7^20),
        b := [ 1 + O(7^20), 1 + O(7^20), 1 + O(7^20), -1 + O(7^20), 2 + O(7^20) 
            ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := -7253842390692000 + O(7^20),
        b := [ 1 + O(7^20), -21761527172076000 + O(7^20), 15826565216055273 + 
            O(7^20), -32312570649446182 + O(7^20), 14507684781384000 + O(7^20) 
        ],
        inf := false>,
    rec<recformat<x, b, inf, xt, bt, index> | 
        x := 0,
        b := [ 1 + O(7^20), 0, 0, 0, 0 ],
        inf := true>
]
*/

for i:=2 to #Qpoints do
 coleman_integrals_on_basis(Qpoints[1],Qpoints[i],data:e:=75);
end for;
/*
(O(7^9) O(7^9) O(7^9) O(7^9) 11452342 + O(7^9) 2808875*7^-1 + O(7^9) 1314185 + 
    O(7^9) 115220003*7^-1 + O(7^9))
9
(O(7^17) O(7^17) O(7^17) O(7^17) 11077637977367*7 + O(7^17) 723737953270644*7^-1
    + O(7^17) -76532500 + O(7^17) 1050159000*7^-1 + O(7^17))
17
(O(7^17) O(7^17) O(7^17) O(7^17) O(7^17) -301560683580935*7^-1 + O(7^17) O(7^17)
    1925291500*7^-1 + O(7^17))
17
(O(7^17) O(7^17) O(7^17) O(7^17) 11077637977367*7 + O(7^17) 603114961058870*7^-1
    + O(7^17) -76532500 + O(7^17) 875132500*7^-1 + O(7^17))
17
(O(7^17) O(7^17) O(7^17) O(7^17) -38772004527972 + O(7^17) 663426457164757*7^-1 
    + O(7^17) -29078694770354 + O(7^17) 962645750*7^-1 + O(7^17))
17
*/
