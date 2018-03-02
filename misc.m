//////////////////////////////////////////////////////////
// Some extra functions not included in Coleman-1.1 yet //  
//////////////////////////////////////////////////////////

eval_R:=function(f,pt,r)

  // Evaluate an element of R at x=x0, y=y0.

  R:=Parent(f); S:=BaseRing(R); Ox:=BaseRing(S); O:=BaseRing(Ox);

  zR:=(Ox!r)/LeadingCoefficient(r);

  x0:=O!pt[1];
  y0:=O!pt[2];
  z0:=Evaluate(zR,x0);
  ev:=O!0;
  C:=Coefficients(f);
  for i:=1 to #C do
    val:=Valuation(C[i]);
    D:=Coefficients(C[i]);
    for j:=1 to #D do
      ev:=ev+Evaluate(D[j],x0)*z0^(val+j-1)*y0^(i-1);
    end for;
  end for;

  return ev;
end function;


eval_mat_R:=function(A,pt,r)

  // Evaluate a matrix over R at x=x0, y=y0.

  R:=BaseRing(A); S:=BaseRing(R); Ox:=BaseRing(S); O:=BaseRing(Ox);

  B:=ZeroMatrix(O,NumberOfRows(A),NumberOfColumns(A));
  for i:=1 to NumberOfRows(A) do
    for j:=1 to NumberOfColumns(A) do
      B[i,j]:=eval_R(A[i,j],pt,r);
    end for;
  end for;

  return B;
end function;


compute_F:=function(Q,W0,Winf,f0,finf,fend)

  // Given functions f0,finf and fend, as vectors of coefficients w.r.t. b^0,b^inf,b^0 respectively, 
  // return f0+finf+fend as a vector w.r.t. b^0 (so convert finf from b^inf to b^0 and take the sum). 
  
  d:=Degree(Q);
  W:=Winf*W0^(-1);

  Qxzzinvd:=Parent(f0);
  Qxzzinv:=BaseRing(Qxzzinvd);
  x1:=Qxzzinv!(BaseRing(Qxzzinv).1);
  Qxxinv:=LaurentSeriesRing(RationalField());

  conv:=Qxzzinvd!Evaluate(Evaluate(finf,Qxxinv.1)*Evaluate(W,Qxxinv.1),x1); // finf converted to basis b^0
  F:=f0+conv+fend;

  return F;
end function;


Qxzzinvd_to_R:=function(f,Q,p,r,R,W0)

  // Convert from Q[x,z,1/z]^d to R, using the basis b^0.

  d:=Degree(Q);

  S:=BaseRing(R);
  Ox:=BaseRing(S);
  y:=R.1;
  z:=S.1;
  x:=Ox.1;

  rQx:=Parent(Numerator(W0[1,1]))!r;
  lc:=LeadingCoefficient(rQx);
  rQx:=rQx/lc; // rQx now corresponds to z=r/LeadingCoefficient(r)

  ordrW0:=ord_r_mat(W0,rQx);

  b0:=[];
  for i:=1 to d do
    b0i:=R!0;
    for j:=1 to d do
      b0i:=b0i+z^(ordrW0)*(R!Ox!(Numerator(W0[i,j]*rQx^(-ordrW0))))*y^(j-1);
    end for;
    b0:=Append(b0,b0i);
  end for;

  f_R:=R!0;  
  for i:=1 to d do
    if f[i] ne 0 then
      for j:=Valuation(f[i]) to Degree(f[i]) do
        coef:=Ox!Coefficient(f[i],j);
        f_R:=f_R+coef*z^j*b0[i]; 
      end for;
    end if;
  end for;

  return f_R;
end function;
