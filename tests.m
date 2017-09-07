load "coleman.m";

print "";
print "//////////////////////////";
print "// 1. An elliptic curve //";
print "//////////////////////////\n";

Q:=y^2-(x^3-10*x+9);
p:=5;
N:=10;
data:=coleman_data(Q,p,N);

print "/////////////////////";
print "// 1a. good points //";
print "/////////////////////\n";

P1:=set_point(0,3,data);  // not bad
P2:=set_point(8,21,data); // not bad
C:=coleman_integrals_on_basis(P1,P2,data);
C;
print "(-117821*5^2 + O(5^10) 1503616 + O(5^10))\n";

print "//////////////////////////";
print "// 1b. finite bad point //";
print "//////////////////////////\n";

P3:=set_bad_point(1,[1,0],false,data); // finite bad
C:=coleman_integrals_on_basis(P1,P3,data:delta:=51);
C;
print "(O(5^10) O(5^10))\n";

print "////////////////////////";
print "// 1c. infinite point //";
print "////////////////////////\n";

P4:=set_bad_point(0,[1,0],true,data); // infinite
C:=coleman_integrals_on_basis(P2,P4,data:delta:=11); 
C;
print "(117821*5^2 + O(5^10) 449509 + O(5^9))\n";


print "///////////////////////////////////////////////";
print "// 2. Bruin-Poonen-Stoll (Winf not diagonal) //";
print "///////////////////////////////////////////////\n";

Q:=y^3 + (-x^2 - 1)*y^2 - x^3*y + x^3 + 2*x^2 + x;
p:=7;
N:=10;
Qp:=pAdicField(p,N);
data:=coleman_data(Q,p,N);

print "/////////////////////";
print "// 2a. good points //";
print "/////////////////////\n";

P1:=set_point(5,-32582624253112412,data);
P2:=set_point(12,25123732258943588,data);
C:=coleman_integrals_on_basis(P1,P2,data);
T:=tiny_integrals_on_basis(P1,P2,data);
C-T;
print "(O(7^11) O(7^10) O(7^10) O(7^10) O(7^9) O(7^9))\n";

print "///////////////////////////";
print "// 2b. finite bad points //";
print "///////////////////////////\n";

P3:=set_point(0,0,data); // finite bad
P4:=set_point(0,1,data); // finite bad
C:=coleman_integrals_on_basis(P1,P3,data:delta:=100); 
C;
print "(-2234847*7 + O(7^9) -15691180 + O(7^9) 17383473 + O(7^9) 651834*7 + O(7^9) 1722235*7^-1 + O(7^8) -6198545*7^-1 + O(7^8))\n";
C:=coleman_integrals_on_basis(P1,P4,data:delta:=100);
C;
print "(712948*7 + O(7^9) -12363515 + O(7^9) -18622514 + O(7^9) 7313571 + O(7^9) -7554227*7^-1 + O(7^8) 18877150*7^-1 + O(7^8))\n";

print "//////////////////////////////////////////";
print "// 2c. infinite points (all unramified) //";
print "//////////////////////////////////////////\n";

P5:=set_bad_point(0,[1,0,1],true,data); // infinite
P6:=set_bad_point(0,[1,1,1],true,data); // infinite
P7:=set_bad_point(0,[1,0,0],true,data); // infinite

C:=coleman_integrals_on_basis(P1,P5,data:delta:=100);
C;
print "(-162150*7^2 + O(7^9) -9516889 + O(7^9) 13528894 + O(7^9) 6856492 + O(7^9) 19063851*7^-1 + O(7^8) 12436440*7^-1 + O(7^8))\n";
C:=coleman_integrals_on_basis(P1,P6,data:delta:=100);
C;
print "(2869652*7 + O(7^9) -2871461 + O(7^9) -8739546 + O(7^9) -14556438 + O(7^9) 2321647*7^-1 + O(7^8) 3349709*7^-1 + O(7^8))\n";
C:=coleman_integrals_on_basis(P1,P7,data:delta:=100);
C;
print "(273453*7 + O(7^9) 8996082 + O(7^9) 11355084 + O(7^9) -16214414 + O(7^9) 14858352*7^-1 + O(7^8) -9306089*7^-1 + O(7^8))\n";


print "/////////////////////////////////////////////////////";
print "// 3. modular curve X0(44) (plane model singular!) //";
print "/////////////////////////////////////////////////////\n";

Q:=y^5+12*x^2*y^3-14*x^2*y^2+(13*x^4+6*x^2)*y-(11*x^6+6*x^4+x^2);
p:=7;
N:=10;
data:=coleman_data(Q,p,N);

print "/////////////////////";
print "// 3a. good points //";
print "/////////////////////\n";

P1:=set_point(1,1,data);
P2:=set_point(8,29647929146699830,data);
C:=coleman_integrals_on_basis(P1,P2,data);
T:=tiny_integrals_on_basis(P1,P2,data);
C-T;
print "(O(7^10) O(7^11) O(7^10) O(7^10) O(7^10) O(7^8) O(7^10) O(7^8))\n";

print "//////////////////////////";
print "// 3b. finite bad point //";
print "//////////////////////////\n";

P3:=set_bad_point(0,[1,0,0,0,0],false,data); // finite bad
C:=coleman_integrals_on_basis(P1,P3,data:delta:=100); // (15(P3-P1)=O on the Jacobian!)
C;
print "(3*7^9 + O(7^10) O(7^9) O(7^9) O(7^9) -10971845 + O(7^9) 4086575*7^-1 + O(7^8) -16736670 + O(7^9) -4312348*7^-1 + O(7^8))\n";

print "////////////////////////"; 
print "// 3c. infinite point //";
print "////////////////////////\n";

P4:=set_bad_point(0,[1,0,0,0,0],true,data); // infinite
C:=coleman_integrals_on_basis(P1,P4,data:delta:=100); // (15(P4-P1)=O on the Jacobian!)
C;
print "(O(7^9) O(7^9) O(7^9) O(7^9) -11056758 + O(7^9) 4086575*7^-1 + O(7^8) -1605787 + O(7^9) -4312348*7^-1 + O(7^8))\n";
