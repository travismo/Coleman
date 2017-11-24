//////////////////////////////////////////////////////
// Example from https://arxiv.org/abs/1007.4617     //
// "On elliptic curves with an isogeny of degree 7" //
//////////////////////////////////////////////////////

load "coleman.m";
Q:=y^7-(x^3-2*x^2-x+1)*(x^3-x^2-2*x+1)^6;
N:=15;

p:=2;

data:=coleman_data(Q,p,N);
Qpoints:=Q_points(data,1000);
L,v:=effective_chabauty(data:Qpoints:=Qpoints,e:=10*p);
print p,L,v;

// Take one point P in each residue disk:

Mlist:=[];

for P in [Qpoints[3],Qpoints[4],Qpoints[6]] do

xt,bt:=local_coord(P,25,data);
Qtt:=LaurentSeriesRing(RationalField());
F2:=FiniteField(2);
F2tt:=LaurentSeriesRing(F2);
d:=Degree(Q);

// Expand each annihilating differential at P modulo 2

omega:=[];
for i:=1 to #v do
  omega[i]:=Qtt!0;
  for j:=1 to 12 do
      for k:=1 to d do 
        omega[i]:=omega[i]+(v[i][j])*Evaluate(data`basis[j][k],Qtt!xt)*Derivative(Qtt!xt)/(-Evaluate(data`r,Qtt!xt));
      end for;
  end for;
  omega[i]:=omega[i]/GCD([IntegerRing()!c : c in Coefficients(omega[i])]);
end for;
omega:=Eltseq(ChangeRing(Vector(omega),F2tt));

// Compute the subspace of the annihilating differentials with order >=2 at P modulo 2

M:=ZeroMatrix(F2,8,2);
for i:=1 to 8 do
  for j:=1 to 2 do
    M[i,j]:=Coefficient(omega[i],j-1);
  end for;
end for;

Mlist:=Append(Mlist,M);

end for;

Dimension(Kernel(Mlist[1]) meet Kernel(Mlist[2]) meet Kernel(Mlist[3]));
// 6

// The dimension is 6 as observed in the paper.

