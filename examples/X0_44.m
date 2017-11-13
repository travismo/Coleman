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
(7940*7^15 + O(7^20) 1985*7^15 + O(7^20) -125*7^15 + O(7^19) 45429*7^14 + 
    O(7^20) -10831335573106609 + O(7^20) -2593400738656141*7^-1 + O(7^18) 
    9242230517248601 + O(7^20) -2924496927479138*7^-1 + O(7^18))
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
(-7^19 + O(7^20) 83*7^16 + O(7^19) O(7^19) 99*7^16 + O(7^19) -5655137081123041 +
    O(7^19) -4221814336566590*7^-1 + O(7^18) 967909028490801 + O(7^19) 
    962645750*7^-1 + O(7^18))
11
(-3*7^19 + O(7^20) -7^20 + O(7^21) -2*7^19 + O(7^20) O(7^20) 3799631722911881*7 
    + O(7^20) -2533089242550254*7^-1 + O(7^18) -76532500 + O(7^20) 
    1050159000*7^-1 + O(7^18))
17
(-16*7^19 + O(7^21) O(7^21) 7^19 + O(7^20) O(7^21) O(7^20) 2955266512239963*7^-1
    + O(7^18) O(7^20) 1925291500*7^-1 + O(7^18))
17
(7^19 + O(7^20) -7^20 + O(7^21) 3*7^19 + O(7^20) O(7^20) 3799631722911881*7 + 
    O(7^20) 5488355754790217*7^-1 + O(7^18) -76532500 + O(7^20) 875132500*7^-1 +
    O(7^18))
17
(O(7^19) O(7^19) O(7^19) O(7^19) -1899816116425628 + O(7^19) 
    -4221814336566590*7^-1 + O(7^18) -1424861778693596 + O(7^19) 962645750*7^-1 
    + O(7^18))
17
*/
