max_prec:=function(Q,p,N,g,W0,Winf,e0,einf);

  // Compute the p-adic precision required for provable correctness

  d:=Degree(Q);
  W:=Winf*W0^(-1);
  
  Nmax:=N+Floor(log(p,-p*(ord_0_mat(W)+1)*einf));
  while (Nmax-Floor(log(p,p*(Nmax-1)*e0))-Floor(log(p,-(ord_inf_mat(W^(-1))+1)*einf)) lt N) do 
    Nmax:=Nmax+1;
  end while;

  Nmax:=Maximum(Nmax,2); 

  return Nmax; // precision to be used in the computations
end function;


frobmatrix:=function(Q,p,N,Nmax,g,r,W0,Winf,G0,Ginf,frobmatb0r,red_list_fin,red_list_inf,basis,integrals,quo_map);

  // Compute the matrix of F_p on H^1(X) mod p^N with respect to 'basis'.

  F:=ZeroMatrix(Rationals(),#basis,#basis);  
  f0list:=[];
  finflist:=[];
  fendlist:=[];

  for i:=1 to #basis do

    dif:=frobenius(basis[i],Q,p,Nmax,r,frobmatb0r);
    dif:=convert_to_Qxzzinvd(dif,Q);

    coefs,f0,finf,fend:=reduce_with_fs(dif,Q,p,N,Nmax,r,W0,Winf,G0,Ginf,red_list_fin,red_list_inf,basis,integrals,quo_map);

    for j:=1 to #basis do
      F[i,j]:=coefs[j];
    end for;
    
    f0list:=Append(f0list,f0);
    finflist:=Append(finflist,finf);
    fendlist:=Append(fendlist,fend);

  end for;
 
  return F,f0list,finflist,fendlist;

end function;


coleman_data:=function(Q,p,N:useU:=false)

  // Takes a polynomial Q in two variables x,y over the rationals which is monic in y.
  // Returns the Coleman data of (the projective nonsingular model of) the curve defined
  // by Q at p to p-adic precision N.

  g:=genus(Q,p);
  r,Delta,s:=auxpolys(Q);

  if not smooth(r,p) then
    error "bad prime";
  end if;

  W0:=mat_W0(Q);
  Winf:=mat_Winf(Q);

  G:=con_mat(Q,Delta,s);
  G0:=W0*Evaluate(G,Parent(W0[1,1]).1)*W0^(-1)+ddx_mat(W0)*W0^(-1);
  Ginf:=Winf*Evaluate(G,Parent(Winf[1,1]).1)*Winf^(-1)+ddx_mat(Winf)*Winf^(-1);

  Jinf,Tinf,Tinfinv:=jordan_inf(Ginf);
  J0,T0,T0inv:=jordan_0(r,G0);
  e0,einf:=ram(J0,Jinf);

  basis,integrals,quo_map:=basis_coho(Q,p,r,W0,Winf,G0,Ginf,J0,Jinf,T0inv,Tinfinv,useU);

  Nmax:=max_prec(Q,p,N,g,W0,Winf,e0,einf);

  frobmatb0r:=froblift(Q,p,Nmax-1,r,Delta,s,W0);

  red_list_fin,red_list_inf:=red_lists(Q,p,Nmax,r,W0,Winf,G0,Ginf,e0,einf,J0,Jinf,T0,Tinf,T0inv,Tinfinv);

  F,f0list,finflist,fendlist:=frobmatrix(Q,p,N,Nmax,g,r,W0,Winf,G0,Ginf,frobmatb0r,red_list_fin,red_list_inf,basis,integrals,quo_map);

  // formatting the output into a record:

  format:=recformat<Q,p,N,W0,Winf,r,Delta,s,G0,Ginf,e0,einf,basis,quo_map,integrals,F,f0list,finflist,fendlist,Nmax,red_list_fin,red_list_inf>;
  out:=rec<format|>;
  out`Q:=Q; out`p:=p; out`N:=N; out`W0:=W0; out`Winf:=Winf; out`r:=r; out`Delta:=Delta; out`s:=s; out`G0:=G0; out`Ginf:=Ginf; 
  out`e0:=e0; out`einf:=einf; out`basis:=basis; out`quo_map:=quo_map; out`integrals:=integrals; out`F:=F; out`f0list:=f0list; 
  out`finflist:=finflist; out`fendlist:=fendlist; out`Nmax:=Nmax; out`red_list_fin:=red_list_fin; out`red_list_inf:=red_list_inf;

  return out;

end function;


set_point:=function(x0,y0,data)

  // Constructs a point from affine coordinates x0,y0. 

  Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0;
  d:=Degree(Q);

  if x0 in RationalField() then
    K:=pAdicField(p,N);
  else
    K:=Parent(x0);
  end if;
  x0:=K!x0; y0:=K!y0;

  if Valuation(x0) lt 0 then
    error "x0 has negative valuation";
  end if;
  
  if (not(W0 eq IdentityMatrix(BaseRing(W0),d))) and (Valuation(Evaluate(data`r,x0)) gt 0) then
    error "W0 is not the identity and r(x0) is zero mod p";
  end if;
  
  format:=recformat<x,b,inf,xt,bt,index>;
  P:=rec<format|>;
  P`inf:=false;
  P`x:=x0;

  y0powers:=[];
  y0powers[1]:=K!1;
  for i:=2 to d do
    y0powers[i]:=(y0)^(i-1);
  end for;
  y0powers:=Vector(y0powers);
  W0x0:=Transpose(Evaluate(W0,x0));

  P`b:=y0powers*W0x0; // the values of the b_i^0 at P

  return P;
end function;


set_bad_point:=function(x,b,inf,data)

  Q:=data`Q; p:=data`p; N:=data`N; 
  Qp:=pAdicField(p,N); d:=Degree(Q);

  format:=recformat<x,b,inf,xt,bt,index>;
  P:=rec<format|>;
  P`inf:=inf;
  P`x:=Qp!x;
  P`b:=[Qp!b[i]:i in [1..d]];

  return P; 

end function;


is_bad:=function(P,data)

  // check whether the point P is bad

  x0:=P`x; r:=data`r;

  if P`inf then // infinite point
    return true;
  end if;

  if Valuation(Evaluate(r,x0)) gt 0 then // finite bad point
    return true;
  end if;

  return false;
  
end function;


is_very_bad:=function(P,data)

  // check whether the point P is very bad

  x0:=P`x; r:=data`r; N:=data`N;

  if P`inf then // infinite point
    if Valuation(x0) ge N then // infinite very bad point
      return true;
    end if;
  else // finite point
    if Valuation(Evaluate(r,x0)) ge N then // finite very bad point
      return true;
    end if;
  end if;

  return false;

end function;


lie_in_same_disk:=function(P1,P2,data)

  // check whether two points P1,P2 lie in the same residue disk

  x1:=P1`x; b1:=P1`b; x2:=P2`x; b2:=P2`b; Q:=data`Q;
  d:=Degree(Q);
  
  if P1`inf ne P2`inf then // one point infinite, other not
    return false;
  else
    for i:=1 to d do
        if Valuation(b1[i]-b2[i]) lt 1 then
          return false;
        end if;
      end for;
      return true;
  end if;
 
end function;


minpoly:=function(f1,f2)

  // computes the minimum polynomial of f2 over Q(f1), where
  // f1,f2 are elements of a 1 dimensional functionfield over Q

  FF:=Parent(f1);
  d1:=Degree(f1);
  d2:=Degree(f2);

  S:=[];
  for i:=0 to d1 do
    for j:=0 to d2 do
      S:=Append(S,f1^j*f2^i);
    end for;
  end for;

  denom:=1;
  for i:=1 to #S do
    E:=Eltseq(S[i]);
    for j:=1 to #E do
      denom:=LCM(denom,Denominator(E[j]));
    end for;
  end for;
  
  maxdeg:=0;
  for i:=1 to #S do
    E:=Eltseq(S[i]);
    for j:=1 to #E do
      deg:=Degree(Numerator(denom*E[j]));
      if  deg gt maxdeg then
        maxdeg:=deg;
      end if;
    end for;
  end for;

  T:=[];
  for i:=1 to #S do
    E:=Eltseq(S[i]);
    v:=[];
    for j:=1 to #E do
      for k:=0 to maxdeg do
        v[(j-1)*(maxdeg+1)+k+1]:=Coefficient(Numerator(denom*E[j]),k);
      end for;  
    end for;
    T:=Append(T,v);
  end for;

  Qx:=PolynomialRing(RationalField());
  Qxy:=PolynomialRing(Qx);

  R:=Basis(NullSpace(Matrix(T)))[1];
  
  poly:=Qxy!0;
  for i:=0 to d1 do
    for j:=0 to d2 do
      poly:=poly+R[i*(d2+1)+j+1]*Qx.1^j*Qxy.1^i;
    end for;
  end for;

  fac:=Factorisation(poly);
  for i:=1 to #fac do
    factor:=fac[i][1];
    test:=FF!0;
    for j:=0 to Degree(factor) do
      test:=test+Evaluate(Coefficient(factor,j),f1)*f2^j;
    end for;
    if test eq 0 then
      poly:=factor;
    end if;
  end for;
  
  return poly;
 
end function;


frobenius_pt:=function(P,data);

  // Computes the image of P under Frobenius

  x0:=P`x; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf;
  d:=Degree(Q); K:=Parent(x0); Ky:=PolynomialRing(K);

  x0p:=x0^p;
  b:=P`b;

  Qx:=RationalFunctionField(RationalField()); Qxy:=PolynomialRing(Qx);

  f:=Qxy!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      f:=f+Coefficient(Coefficient(Q,i),j)*Qxy.1^i*Qx.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  if not is_bad(P,data) then // finite good point
    
    W0invx0:=Transpose(Evaluate(W0^(-1),x0));

    ypowers:=Vector(b)*W0invx0;
    y0:=ypowers[2];
  
    C:=Coefficients(Q);
    D:=[];
    for i:=1 to #C do
      D[i]:=Evaluate(C[i],x0p);
    end for;
    fy:=Ky!D;

    y0p:=HenselLift(fy,y0^p); // Hensel lifting
  
    y0ppowers:=[];
    y0ppowers[1]:=K!1;
    for i:=2 to d do
      y0ppowers[i]:=y0p^(i-1);
    end for;
    y0ppowers:=Vector(y0ppowers);

    W0x0:=Transpose(Evaluate(W0,x0));
  
    b:=y0ppowers*W0x0;

  elif P`inf then // infinite point
  
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+Winf[i,j]*FF.1^(j-1);
      end for;
      poly:=minpoly(FF!(1/Qx.1),bi);

      C:=Coefficients(poly);
      D:=[];
      for i:=1 to #C do
        D[i]:=Evaluate(C[i],x0p); 
      end for;
      fy:=Ky!D;

      zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
      done:=false;
      j:=1;
      while not done and j le #zeros do
        if Valuation(zeros[j][1]-b[i]^p) gt p then // was gt 0 before 
          done:=true;
          b[i]:=zeros[j][1];
        end if;
        j:=j+1;
      end while;
      if not done then
        error "Frobenius does not converge at P";
      end if;
    end for;

  else // finite bad point

   for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+W0[i,j]*FF.1^(j-1);
      end for;
      poly:=minpoly(FF!Qx.1,bi);

      C:=Coefficients(poly);
      D:=[];
      for i:=1 to #C do
        D[i]:=Evaluate(C[i],x0p); 
      end for;
      fy:=Ky!D;

      zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
      done:=false;
      j:=1;
      while not done and j le #zeros do
        if Valuation(zeros[j][1]-b[i]^p) gt p then // was gt 0 before
          done:=true;
          b[i]:=zeros[j][1];
        end if;
        j:=j+1;
      end while;
      if not done then
        error "Frobenius does not converge at P";
      end if;
    end for;

  end if;
  
    P`x:=x0p;
    P`b:=b;

  return P;
end function;


local_data:=function(P,data)

  // For a point P, returns the ramification index of the map x on the residue disk at P

  Q:=data`Q; p:=data`p; W0:=data`W0; Winf:=data`Winf; x0:=P`x; b:=P`b; d:=Degree(Q);

  if not is_bad(P,data) then
    eP:=1;
  else     
    Fp:=FiniteField(p); Fpx:=RationalFunctionField(Fp); Fpxy:=PolynomialRing(Fpx);
    f:=Fpxy!0;
    for i:=0 to d do
      for j:=0 to Degree(Coefficient(Q,i)) do
        f:=f+(Fp!Coefficient(Coefficient(Q,i),j))*Fpxy.1^i*Fpx.1^j;
      end for;
    end for;  
    FFp:=FunctionField(f); // function field of curve mod p
    
    if P`inf then
      places:=InfinitePlaces(FFp); // infinite places of function field of curve mod p
      W:=Winf;
    else
      Px0:=Zeros(Fpx.1-Fp!x0)[1]; 
      places:=Decomposition(FFp,Px0); // places of function field of curve mod p lying over x0 mod p
      W:=W0;
    end if;

    bmodp:=[]; // elements of b mod p, where b is either b^0 or b^inf
    for i:=1 to d do
      f:=FFp!0;
      for j:=1 to d do
        f:=f+(Fpx!W[i,j])*FFp.1^(j-1);
      end for;
      bmodp[i]:=f;
    end for;

    // compute values of b_j at places 
      
    values:=[];
    for i:=1 to #places do
      placei:=places[i];
      bjplacei:=[];
      for j:=1 to d do
        bjplacei[j]:=Evaluate(bmodp[j],placei); // value of b^_j at i-th place 
      end for;
      values:=Append(values,bjplacei);
    end for;

    // identify the right place

    v:=[];
    for j:=1 to d do
      v[j]:=Fp!b[j];
    end for;
    
    done:=false;
    i:=1;
    while (not done) and (i le #places) do  
      if v eq values[i] then // P reduces to places[i] mod p
        place:=places[i];
        done:=true;
      end if;
      i:=i+1;
    end while;

    if not done then
      error "Point does not lie on curve";
    end if;

    eP:=RamificationIndex(place);

    if eP eq 1 then
      index:=0;
    else
      done:=false;
      i:=1;
      while not done do
        ord:=Valuation(bmodp[i]-Evaluate(bmodp[i],place),place);
        if ord eq 1 then
          index:=i;
          done:=true;
        end if;
        i:=i+1;
      end while;
    end if;
  
  end if;

  return eP,index,place,bmodp;

end function;


hensel_lift:=function(fy,root);

  // Finds a root of the polynomial fy over Qp[[t]]
  // by Hensel lifting from an approximate root.
  //
  // Assumes that the starting criterion for Hensel's 
  // lemma is satisfied

  Kty:=Parent(fy);
  Kt:=BaseRing(Kty);
  K:=BaseRing(Kt);
  tprec:=Precision(Kt); // t-adic precision
  Qt:=PowerSeriesRing(RationalField(),tprec);
  Qty:=PolynomialRing(Qt);
  p:=Prime(K);
  pprec:=Precision(K);  // p-adic precision
  Zp:=pAdicRing(p,pprec);
  Zpt:=PowerSeriesRing(Zp,tprec);  

  fy:=Qty!fy;
  derfy:=Derivative(fy);  

  v1:=Valuation(Qt!Zpt!Evaluate(fy,root));
  v2:=Valuation(Qt!Zpt!Evaluate(derfy,root));

  if not v1 gt 2*v2 then
    error "Condition Hensel's Lemma not satisfied";
  end if;

  if Qt!root eq 0 then
    deg_root:=1;
  else
    deg_root:=Degree(Qt!root);
  end if;

  prec_seq:=[];
  k:=tprec;
  
  while k gt v1 do
    prec_seq:=Append(prec_seq,k);
    k:=Ceiling(k/2+v2);
  end while;
  prec_seq:=Reverse(prec_seq);

  for j:=1 to #prec_seq do
    root:=Qt!root;
    root:=ChangePrecision(root,prec_seq[j]);
    root:=root-Evaluate(fy,root)/Evaluate(derfy,root);
    root:=Kt!root;
  end for;

  return root;

end function;


mod_p_prec:=function(fy)

  // Finds the t-adic precision necessary to separate the roots
  // of the polynomial fy over Qp[[t]] modulo p.

  Kty:=Parent(fy);
  Kt:=BaseRing(Kty);
  tprec:=Precision(Kt);
  K:=BaseRing(Kt);
  p:=Prime(K);
  Fp:=FiniteField(p);
  Fpt:=PowerSeriesRing(Fp,tprec);
  Fpty:=PolynomialRing(Fpt);

  fymodp:=Fpty!fy;
  derfymodp:=Derivative(fymodp);

  done:=false;
  prec:=1;

  while not done do
    prec:=prec+1;
    zeros:=Roots(fymodp,prec); // can run into trouble if prec is too high compared to tprec...
    done:=true;
    for i:=1 to #zeros do
      root:=zeros[i][1];
      if Evaluate(fymodp,root) eq 0 then
        v1:=prec;
      else
        v1:=Valuation(Evaluate(fymodp,root));
      end if;
      v2:=Valuation(Evaluate(derfymodp,root));
      if zeros[i][2] gt 1 or v1 le 2*v2 then
        done:=false;
      end if;
    end for;
  end while;

  modpprec:=Maximum([Valuation(z[1]):z in zeros])+prec;

  return modpprec;

end function;


approx_root:=function(fy,y0,modpprec,expamodp)

  // Computes an approximation to t-adic precision modpprec of 
  // a root of the polynomial fy over Qp[[t]] which is congruent to:
  //
  // y0 modulo t
  // expamodp modulo p 
  //
  // This approximation is then used as root in hensel_lift.

  Kty:=Parent(fy);
  Kt:=BaseRing(Kty);
  tprec:=Precision(Kt); // t-adic precision
  K:=BaseRing(Kt);

  p:=Prime(K);
  Fp:=FiniteField(p);
  pprec:=Precision(K);  // p-adic precision
  Zp:=pAdicRing(p,pprec);
  Zpt:=PowerSeriesRing(Zp,tprec);
  Zpz:=PolynomialRing(Zp);

  Qt:=PowerSeriesRing(RationalField(),tprec);
  Qty:=PolynomialRing(Qt);
  Qz:=PolynomialRing(RationalField());
  Qzt:=PowerSeriesRing(Qz,tprec);
  
  root:=Kt!y0;
  for i:=1 to modpprec-1 do
    newroot:=root+Kty.1*Kt.1^i;
    C:=Coefficients(fy);
    fynewroot:=Kty!0;
    for i:=1 to #C do
       fynewroot:=fynewroot+(Kt!C[i])*newroot^(i-1);
    end for;
    fynewroot:=Qty!Kty!fynewroot;
    fznewroot:=Qzt!0;
    for j:=0 to Degree(fynewroot) do
      for k:=0 to tprec-1 do
        fznewroot:=fznewroot+Coefficient(Coefficient(fynewroot,j),k)*(Qz.1)^j*(Qzt.1)^k;
      end for;
    end for;
    zeros:=Roots(Zpz!Coefficient(fznewroot,Valuation(fznewroot)));
    count:=0;
    for j:=1 to #zeros do
      sol:=zeros[j][1];
      if Fp!sol eq Coefficient(expamodp,i) then
        count:=count+1;
        coef:=sol;
      end if;
    end for;
    if count ne 1 then // I think this should never happen, but let's keep this for now.
      error "Coefficients in expansion of approximation coincide in Fp but not in Qp";
    end if;
    root:=Evaluate(newroot,coef);
  end for;
  
  return root;

end function;


mod_p_expansion:=function(f,place,tmodp,modpprec);

  // Finds the power series expansion of f in the function field
  // modulo p at place with respect to local parameter tmodp to
  // absolute precision modpprec.

  FFp:=Parent(f);
  Fpx:=BaseRing(FFp);
  Fp:=BaseRing(Fpx);
  Fpt:=PowerSeriesRing(Fp,modpprec);

  expamodp:=Fpt!0;
  for i:=0 to modpprec-1 do
    val:=Evaluate(f,place);
    expamodp:=expamodp+val*Fpt.1^i;
    f:=(f-val)/tmodp;
  end for;
  
  return expamodp;
  
end function;


local_coord:=function(P,prec,data);

  // Computes powerseries expansions of x and
  // the b^0_i or b^infty_i (depending on whether
  // P is infinite or not) in terms of the local
  // coordinate computed by local_data.

  if assigned P`xt and Precision(Parent(P`xt)) ge prec then
    xt:=P`xt;
    bt:=P`bt;
    index:=P`index;
    return xt,bt,index;
  end if;

  x0:=P`x; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; d:=Degree(Q); b:=P`b;
  K:=Parent(x0); Kt<t>:=PowerSeriesRing(K,prec); Kty:=PolynomialRing(Kt);
  Qx:=RationalFunctionField(RationalField()); Qxy:=PolynomialRing(Qx);
  Fp:=FiniteField(p);

  f:=Qxy!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      f:=f+Coefficient(Coefficient(Q,i),j)*Qxy.1^i*Qx.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  if not is_bad(P,data) then // finite good point

    xt:=t+x0;

    W0invx0:=Transpose(Evaluate(W0^(-1),x0));
    ypowers:=Vector(b)*W0invx0;
    y0:=ypowers[2];

    C:=Coefficients(Q);
    D:=[];
    for i:=1 to #C do
      D[i]:=Evaluate(C[i],xt); 
    end for;
    fy:=Kty!D;
    derfy:=Derivative(fy);

    yt:=hensel_lift(fy,Kt!y0);

    ypowerst:=[];
    ypowerst[1]:=FieldOfFractions(Kt)!1;
    ypowerst[2]:=yt;
    for i:=3 to d do
      ypowerst[i]:=ypowerst[i-1]*yt;
    end for;
    bt:=Eltseq(Vector(ypowerst)*Transpose(Evaluate(W0,xt)));

    btnew:=[];
    for i:=1 to d do
      btnew[i]:=Kt!bt[i];
    end for;
    bt:=btnew;

    index:=0;

  elif P`inf then // infinite point

    eP,index,place,bmodp:=local_data(P,data);
    FFp:=Parent(bmodp[1]);
    Fpx:=BaseRing(FFp);

    bfun:=[];
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+Winf[i,j]*FF.1^(j-1);
      end for;
      bfun[i]:=bi;
    end for;
    
    if eP eq 1 then // P is an infinite point that is not ramified
      
      xt:=t+x0;
      bt:=[];

      for i:=1 to d do

        poly:=minpoly(FF!(1/Qx.1),bfun[i]);

        C:=Coefficients(poly);
        D:=[];
        for j:=1 to #C do
          D[j]:=Evaluate(C[j],xt); 
        end for;
        fy:=Kty!D;
        derfy:=Derivative(fy);

        modpprec:=mod_p_prec(fy);

        if assigned P`bt and Precision(Parent(P`bt[i])) ge modpprec then
          approxroot:=P`bt[i];
        else
          tmodp:=1/Fpx.1-Fp!x0;
          expamodp:=mod_p_expansion(bmodp[i],place,tmodp,modpprec);
          approxroot:=approx_root(fy,b[i],modpprec,expamodp);
        end if;

        bti:=hensel_lift(fy,approxroot);
        bt[i]:=bti;

      end for;

    else // P is an infinite point that is ramified

      poly:=minpoly(bfun[index],FF!1/(Qx.1));

      C:=Coefficients(poly);
      D:=[];
      for j:=1 to #C do
        D[j]:=Evaluate(C[j],t+b[index]); 
      end for;
      fy:=Kty!D;
      derfy:=Derivative(fy);

      modpprec:=mod_p_prec(fy);

      if assigned P`xt and Precision(Parent(P`xt)) ge modpprec then
        approxroot:=P`xt;
      else
        tmodp:=bmodp[index]-Fp!b[index];
        expamodp:=mod_p_expansion(FFp!1/Fpx.1,place,tmodp,modpprec);
        approxroot:=approx_root(fy,x0,modpprec,expamodp);
      end if;

      xt:=hensel_lift(fy,approxroot);

      bt:=[];
      for i:=1 to d do 
      
        if i eq index then
          bt[i]:=t+b[i];
        else
          
          poly:=minpoly(bfun[index],bfun[i]);

          C:=Coefficients(poly);
          D:=[];
          for j:=1 to #C do
            D[j]:=Evaluate(C[j],t+b[index]); 
          end for;

          fy:=Kty!D;
          derfy:=Derivative(fy);

          modpprec:=mod_p_prec(fy);

          if assigned P`bt and Precision(Parent(P`bt[i])) ge modpprec then
            approxroot:=P`bt[i];
          else
            tmodp:=bmodp[index]-Fp!b[index];
            expamodp:=mod_p_expansion(bmodp[i],place,tmodp,modpprec);
            approxroot:=approx_root(fy,b[i],modpprec,expamodp);
          end if;

          bti:=hensel_lift(fy,approxroot);
          bt[i]:=bti;

        end if;
 
      end for;

    end if;

  else // finite bad point

    eP,index,place,bmodp:=local_data(P,data);
    FFp:=Parent(bmodp[1]);
    Fpx:=BaseRing(FFp);

    bfun:=[];
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+W0[i,j]*FF.1^(j-1);
      end for;
      bfun[i]:=bi;
    end for;

    if eP eq 1 then // P is a finite point that is not ramified

      xt:=t+x0;
      bt:=[];
      for i:=1 to d do
        
        poly:=minpoly(FF!Qx.1,bfun[i]);

        C:=Coefficients(poly);
        D:=[];
        for j:=1 to #C do
          D[j]:=Evaluate(C[j],xt); 
        end for;
        fy:=Kty!D;
        derfy:=Derivative(fy);

        modpprec:=mod_p_prec(fy);

        if assigned P`bt and Precision(Parent(P`bt[i])) ge modpprec then
          approxroot:=P`bt[i];
        else
          tmodp:=Fpx.1-Fp!x0;
          expamodp:=mod_p_expansion(bmodp[i],place,tmodp,modpprec);
          approxroot:=approx_root(fy,b[i],modpprec,expamodp);
        end if;

        bti:=hensel_lift(fy,approxroot);
        bt[i]:=bti;

      end for;

    else // P is a finite point that ramifies

      poly:=minpoly(bfun[index],FF!Qx.1);

      C:=Coefficients(poly);
      D:=[];
      for j:=1 to #C do
        D[j]:=Evaluate(C[j],t+b[index]); 
      end for;
      fy:=Kty!D;
      derfy:=Derivative(fy);

      modpprec:=mod_p_prec(fy);

      if assigned P`xt and Precision(Parent(P`xt)) ge modpprec then
        approxroot:=P`xt;
      else
        tmodp:=bmodp[index]-Fp!b[index];
        expamodp:=mod_p_expansion(FFp!Fpx.1,place,tmodp,modpprec);
        approxroot:=approx_root(fy,x0,modpprec,expamodp);
      end if;

      xt:=hensel_lift(fy,Kt!x0); 

      bt:=[];
      for i:=1 to d do 
      
        if i eq index then
          bt[i]:=t+b[i];
        else
          
          poly:=minpoly(bfun[index],bfun[i]);

          C:=Coefficients(poly);
          D:=[];
          for j:=1 to #C do
            D[j]:=Evaluate(C[j],t+b[index]);
          end for;

          fy:=Kty!D;
          derfy:=Derivative(fy);

          modpprec:=mod_p_prec(fy);

          if assigned P`bt and Precision(Parent(P`bt[i])) ge modpprec then
            approxroot:=P`bt[i];
          else
            tmodp:=bmodp[index]-Fp!b[index];
            expamodp:=mod_p_expansion(bmodp[i],place,tmodp,modpprec);
            approxroot:=approx_root(fy,b[i],modpprec,expamodp);
          end if;

          bti:=hensel_lift(fy,approxroot);
          bt[i]:=bti;

        end if;
 
      end for;

    end if;

  end if;

  return xt,bt,index;

end function;


tiny_integral_prec:=function(prec,delta,maxpoleorder,maxdegree,mindegree,val,data);

  // Determines the p-adic precision to which tiny_integrals_on_basis is correct.

  N:=data`N; p:=data`p;

  // Precision loss from terms of positive order we do consider:

  m1:=N*delta-val;
  for i:=1 to maxdegree do
    m1:=Minimum(m1,N*delta+i-delta*Floor(Log(p,i+1)));
  end for;  

  // Precision loss from terms we omit:

  m2:=mindegree+2-delta*Floor(Log(p,mindegree+2));
  for i:=mindegree+2 to Ceiling(delta/Log(p)) do
    m2:=Minimum(m2,i+1-delta*Floor(Log(p,i+1)));
  end for;

  // Precision loss from terms of negative order

  m3:=N*delta-val;
  if maxpoleorder ge 2 then
    m3:=N*delta-val-maxpoleorder*val-delta*Floor(Log(p,maxpoleorder-1));
  end if;

  m:=Minimum([m1,m2,m3]);

  return m/delta;

end function;


find_bad_point_in_disk:=function(P,data);

  // Find the bad point in the residue disk at P.

  x0:=P`x; b:=P`b; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; r:=data`r;
  d:=Degree(Q); K:=Parent(x0); Ky:=PolynomialRing(K);

  if not is_bad(P,data) then
    error "Residue disk does not contain a bad point";
  end if;

  Qx:=RationalFunctionField(RationalField()); Qxy:=PolynomialRing(Qx);

  f:=Qxy!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      f:=f+Coefficient(Coefficient(Q,i),j)*Qxy.1^i*Qx.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  format:=recformat<x,b,inf,xt,bt,index>;
  Pbad:=rec<format|>;
  
  if P`inf then // infinite point
    Pbad`inf:=true;
    x0:=K!0;
    Pbad`x:=x0;
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+Winf[i,j]*FF.1^(j-1);
      end for;
      poly:=minpoly(FF!(1/Qx.1),bi);

      C:=Coefficients(poly);
      D:=[];
      for i:=1 to #C do
        D[i]:=Evaluate(C[i],x0); 
      end for;
      fy:=Ky!D;

      zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
      done:=false;
      j:=1;
      while not done and j le #zeros do
        if Valuation(zeros[j][1]-b[i]) gt 0 then
          done:=true;
          b[i]:=zeros[j][1];
        end if;
        j:=j+1;
      end while;
    end for;
    Pbad`b:=b;
  else // finite bad point
    Pbad`inf:=false;
    rQp:=Ky!Coefficients(r);
    x0:=HenselLift(rQp,x0);
    Pbad`x:=x0;
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+W0[i,j]*FF.1^(j-1);
      end for;
      poly:=minpoly(FF!Qx.1,bi);

      C:=Coefficients(poly);
      D:=[];
      for i:=1 to #C do
        D[i]:=Evaluate(C[i],x0); 
      end for;
      fy:=Ky!D;

      zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
      done:=false;
      j:=1;
      while not done and j le #zeros do
        if Valuation(zeros[j][1]-b[i]) gt 0 then
          done:=true;
          b[i]:=zeros[j][1];
        end if;
        j:=j+1;
      end while;
      Pbad`b:=b;
    end for;
    
  end if;

  return Pbad;

end function;


tiny_integrals_on_basis:=function(P1,P2,data:prec:=0,P:=0);

  // Compute tiny integrals of basis elements from P1 to P2.
  // If P1 is not defined over Qp (but a totally ramified 
  // extension) then a point P defined over Qp in the same
  // residue disk as P1 has to be specified.

  x1:=P1`x; x2:=P2`x; b1:=P1`b; b2:=P2`b; Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; r:=data`r; basis:=data`basis; N:=data`N;
  d:=Degree(Q); lc_r:=LeadingCoefficient(r); W:=Winf*W0^(-1); K:=Parent(x1);

  if not lie_in_same_disk(P1,P2,data) then
    error "the points do not lie in the same residue disk";
  end if;

  if (x1 eq x2) and (b1 eq b2) then
    return RSpace(K,#basis)!0,N*Degree(K);
  end if;

  if Degree(K) gt 1 then // P1 needs to be defined over Qp
    tinyPtoP2,NtinyPtoP2:=$$(P,P2,data);
    tinyPtoP1,NtinyPtoP1:=$$(P,P1,data);
    return tinyPtoP2-tinyPtoP1,Minimum(NtinyPtoP2,NtinyPtoP1);
  end if;

  if not Type(P) eq Rec then
    P:=P1;
  end if;

  if is_bad(P,data) and not is_very_bad(P,data) then // on a bad disk P1 needs to be very bad
    P:=find_bad_point_in_disk(P,data);
    tinyPtoP2,NtinyPtoP2:=$$(P,P2,data);
    tinyPtoP1,NtinyPtoP1:=$$(P,P1,data);
    return tinyPtoP2-tinyPtoP1,Minimum(NtinyPtoP2,NtinyPtoP1);
  end if;

  delta:=Degree(Parent(x2));

  if prec eq 0 then // no t-adic precision specified
    prec:=1;
    while Floor(prec/delta)+1-Floor(Log(p,prec+1)) lt N do
      prec:=prec+1;
    end while;
  end if;
  prec:=Maximum([prec,-2*ord_0_mat(W)*data`einf,50]); // temporary, look at this again

  Kt:=LaurentSeriesRing(K,prec);
  OK:=RingOfIntegers(K);
  OKt:=LaurentSeriesRing(OK,prec);

  xt,bt,index:=local_coord(P1,prec,data);
  
  Qt<t>:=LaurentSeriesRing(RationalField(),prec);
  xt:=Qt!xt;
  btnew:=[Qt|];
  for i:=1 to d do
    btnew[i]:=Qt!bt[i];
  end for;
  bt:=Vector(btnew);

  if P1`inf then
    xt:=1/xt;
    xt:=Qt!Kt!xt; 
    Winv:=W0*Winf^(-1);          
    bt:=bt*Transpose(Evaluate(Winv,xt));
    for i:=1 to d do
      bt[i]:=Qt!(Kt!bt[i]);
    end for; 
  end if;

  if P1`inf or not is_bad(P1,data) then 
    denom:=Qt!Kt!(1/Evaluate(r,xt));
  else
    Qp:=pAdicField(p,N);
    Qpx:=PolynomialRing(Qp);
    rQp:=Qpx!r;
    zero:=HenselLift(rQp,x1);
    sQp:=rQp div (Qpx.1-zero);
    denom:=Qt!Kt!((Qt!OKt!(xt-Coefficient(xt,0)))^(-1)*(Qt!Kt!(1/Evaluate(sQp,xt))));
  end if;

  derxt:=Qt!Kt!Derivative(xt); 
  diffs:=[];
  for i:=1 to #basis do
    basisxt:=Evaluate(basis[i],xt);
    for j:=1 to d do
      basisxt[1][j]:=Qt!Kt!basisxt[1][j];
    end for;
    diffs[i]:=InnerProduct(Vector(basisxt*derxt*lc_r*denom),bt);
    diffs[i]:=Qt!Kt!diffs[i];
    if Coefficient(diffs[i],-1) ne 0 then
      diffs[i]:=diffs[i]-Coefficient(diffs[i],-1)*Kt.1^(-1); // temporary, deal with logs later, important for double integrals
    end if;
  end for;

  maxpoleorder:=-(Minimum([Valuation(diffs[i]): i in [1..#basis]])); 
  maxdegree:=Maximum([Degree(diffs[i]): i in [1..#basis]]);
  mindegree:=Minimum([Degree(diffs[i]): i in [1..#basis]]);

  indefints:=[];
  for i:=1 to #basis do
    indefints := Append(indefints, Integral(diffs[i]));
  end for;

  tinyP1toP2:=[];
  for i:=1 to #basis do
    if index eq 0 then // x-x(P1) is the local coordinate
      tinyP1toP2[i]:=Evaluate(indefints[i],x2-x1);
      val:=Valuation(x2-x1);
    else // b[index]-b[index](P1) is the local coordinate
      tinyP1toP2[i]:=Evaluate(indefints[i],b2[index]-b1[index]);
      val:=Valuation(b2[index]-b1[index]);
    end if;
  end for;

  NtinyP1toP2:=tiny_integral_prec(prec,delta,maxpoleorder,maxdegree,mindegree,val,data);

  return Vector(tinyP1toP2),NtinyP1toP2;

end function;


tiny_integrals_on_basis_to_z:=function(P,data:prec:=0);

  // Compute tiny integrals of basis elements from P to an
  // arbitrary point in its residue disk as a power series
  // in the local parameter there. The series expansions xt
  // and bt of the coordinates on the curve in terms of this 
  // local parameter are also returned.

  x0:=P`x; b:=P`b; Q:=data`Q; p:=data`p; N:=data`N; basis:=data`basis; r:=data`r; W0:=data`W0; Winf:=data`Winf;
  d:=Degree(Q); lc_r:=LeadingCoefficient(r); W:=Winf*W0^(-1); K:=Parent(x0);

  if is_bad(P,data) and not is_very_bad(P,data) then // on a bad disk P needs to be very bad
    P1:=find_bad_point_in_disk(P,data);  
  else
    P1:=P;
  end if;
  x1:=P1`x;

  IPP1,NIPP1:=tiny_integrals_on_basis(P,P1,data:prec:=prec);

  if prec eq 0 then // no precision specified
    prec:=1;
    while Floor(prec)+1-Floor(Log(p,prec+1)) lt N do
      prec:=prec+1;
    end while;
  end if;
  prec:=Maximum([prec,-2*ord_0_mat(W)*data`einf,50]); // temporary, TODO look at this again

  Kt<t>:=LaurentSeriesRing(K,prec);
  OK:=RingOfIntegers(K);
  OKt:=LaurentSeriesRing(OK,prec);

  xt,bt,index:=local_coord(P1,prec,data);

  xtold:=xt;
  btold:=bt;

  Qt<t>:=LaurentSeriesRing(RationalField(),prec);
  xt:=Qt!xt;
  btnew:=[Qt|];
  for i:=1 to d do
    btnew[i]:=Qt!bt[i];
  end for;
  bt:=Vector(btnew);

  if P1`inf then
    xt:=1/xt;
    xt:=Qt!Kt!xt; 
    Winv:=W0*Winf^(-1);          
    bt:=bt*Transpose(Evaluate(Winv,xt));
    for i:=1 to d do
      bt[i]:=Qt!(Kt!bt[i]);
    end for; 
  end if;

  if P1`inf or not is_bad(P1,data) then 
    denom:=Qt!Kt!(1/Evaluate(r,xt));
  else
    Qp:=pAdicField(p,N);
    Qpx:=PolynomialRing(Qp);
    rQp:=Qpx!r;
    zero:=HenselLift(rQp,x1);
    sQp:=rQp div (Qpx.1-zero);
    denom:=Qt!Kt!((Qt!OKt!(xt-Coefficient(xt,0)))^(-1)*(Qt!Kt!(1/Evaluate(sQp,xt))));
  end if;

  derxt:=Qt!Kt!Derivative(xt); 
  diffs:=[];
  for i:=1 to #basis do
    basisxt:=Evaluate(basis[i],xt);
    for j:=1 to d do
      basisxt[1][j]:=Qt!Kt!basisxt[1][j];
    end for;
    diffs[i]:=InnerProduct(Vector(basisxt*derxt*lc_r*denom),bt);
    diffs[i]:=Qt!Kt!diffs[i];
    if Coefficient(diffs[i],-1) ne 0 then
      diffs[i]:=diffs[i]-Coefficient(diffs[i],-1)*t^(-1); // temporary, TODO deal with logs later, important for double integrals
    end if;
  end for;

  indefints:=[];
  for i:=1 to #basis do
    indefints := Append(indefints, Integral(diffs[i]));
  end for;

  xt:=xtold;
  bt:=Vector(btold);

  return Vector(indefints)+IPP1,xt,bt,NIPP1;

end function;


evalf0:=function(f0,P,data);

  // TODO comment
 
  x0:=P`x; b:=P`b; Q:=data`Q; r:=data`r; W0:=data`W0; Winf:=data`Winf; N:=data`N; Nmax:=data`Nmax; p:=data`p;
  d:=Degree(Q); lcr:=LeadingCoefficient(r); K:=Parent(x0);

  Nf0P:=N*Degree(K);

  if P`inf then 
    Winv:=W0*Winf^(-1); 
    Qx:=BaseRing(Winv);
    b:=Vector(b)*Transpose(Evaluate(Evaluate(Winv,1/Qx.1),x0)); // values of the b_i^0 at P
    z0:=Evaluate(r,1/x0)/lcr;
    f0P:=K!0;
    for i:=1 to d do
      f0i:=f0[i];
      C:=Coefficients(f0i);
      val:=Valuation(f0i);
      for j:=1 to #C do
        D:=Coefficients(C[j]);
        for k:=1 to #D do
          f0P:=f0P+(K!D[k])*(1/x0)^(k-1)*z0^(val+j-1)*b[i];
        end for;
      end for;
    end for;
    Nf0P:=N*Degree(K)+(ord_inf_mat(Winv)+1)*Valuation(x0);
  else
    z0:=Evaluate(r,x0)/lcr;  
    f0P:=K!0;
    for i:=1 to d do
      f0i:=f0[i];
      C:=Coefficients(f0i);
      val:=Valuation(f0i);
      for j:=1 to #C do
        D:=Coefficients(C[j]);
        for k:=1 to #D do
          f0P:=f0P+(K!D[k])*x0^(k-1)*z0^(val+j-1)*b[i];
        end for;
      end for;
    end for;
    Nf0P:=N*Degree(K)-p*(Nmax-1)*Valuation(z0); // TODO this is error of terms we did consider, take error of terms we ignored into account as well
  end if;

  return f0P,Nf0P/Degree(K);

end function;


evalfinf:=function(finf,P,data);

  // TODO comment

  x0:=P`x; b:=P`b; Q:=data`Q; W0:=data`W0; Winf:=data`Winf; N:=data`N; p:=data`p;
  d:=Degree(Q); K:=Parent(x0); 

  W:=Winf*W0^(-1); 

  if P`inf then
    finfP:=K!0;
    for i:=1 to d do
      finfi:=finf[i];
      C:=Coefficients(finfi);
      val:=Valuation(finfi);
      for j:=1 to #C do
        finfP:=finfP+(K!C[j])*(1/x0)^(val+j-1)*b[i];
      end for;
    end for;
    NfinfP:=N*Degree(K)+p*(ord_0_mat(W)+1)*Valuation(x0);
  else 
    finf:=finf*ChangeRing(W,BaseRing(finf));
    finfP:=K!0;
    for i:=1 to d do
      finfi:=finf[i];
      C:=Coefficients(finfi);
      val:=Valuation(finfi);
      for j:=1 to #C do
        finfP:=finfP+(K!C[j])*x0^(val+j-1)*b[i];
      end for;
    end for;
    NfinfP:=N*Degree(K);
  end if;

  return finfP, NfinfP/Degree(K);

end function;


evalfend:=function(fend,P,data);

  // TODO comment

  x0:=P`x; b:=P`b; Q:=data`Q; W0:=data`W0; Winf:=data`Winf; N:=data`N;
  d:=Degree(Q);
  K:=Parent(x0);


  if P`inf then
    Winv:=W0*Winf^(-1);
    Qx:=BaseRing(Winv);
    b:=Vector(b)*Transpose(Evaluate(Evaluate(Winv,1/Qx.1),x0)); // values of the b_i^0 at P
    fendP:=K!0;
    for i:=1 to d do
      fendi:=fend[i];
      C:=Coefficients(fendi);
      for j:=1 to #C do
        fendP:=fendP+(K!C[j])*(1/x0)^(j-1)*b[i];
      end for;
    end for;
    NfendP:=N*Degree(K)+(ord_0_mat(Winf)+1)*Valuation(x0);
  else
    fendP:=K!0;
    for i:=1 to d do
      fendi:=fend[i];
      C:=Coefficients(fendi);
      for j:=1 to #C do
        fendP:=fendP+(K!C[j])*x0^(j-1)*b[i];
      end for;
    end for;
    NfendP:=N*Degree(K);
  end if;

  return fendP, NfendP/Degree(K);

end function;


round_to_Qp:=function(L)

  // rounds a vector over a totally ramified extension of Qp to one over Qp

  K:=CoefficientRing(L);
  deg:=Degree(K);
  e:=Precision(K);

  l:=[];
  for i:=1 to #Eltseq(L) do
    l[i]:=Eltseq(L[i])[1];  
    e:=Minimum(e,Valuation(L[i]-l[i]));
  end for;

  return Vector(l),e/deg;

end function;


coleman_integrals_on_basis:=function(P1,P2,data:delta:=1)

  // Integrals of basis elements from P1 to P2. 

  F:=data`F; Q:=data`Q; basis:=data`basis; x1:=P1`x; f0list:=data`f0list; finflist:=data`finflist; fendlist:=data`fendlist; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf;
  d:=Degree(Q); K:=Parent(x1); W:=Winf*W0^(-1);

  prec:=1;
  while Floor(prec/delta)+1-Floor(Log(p,prec+1)) lt N do
    prec:=prec+1;
  end while;
  prec:=Maximum([prec,-2*ord_0_mat(W)*data`einf,50]); // temporary, TODO look at this again

  if is_bad(P1,data) then
    xt,bt,index:=local_coord(P1,prec,data);
    P1`xt:=xt;       // to avoid recomputing
    P1`bt:=bt;       // to avoid recomputing
    P1`index:=index; // to avoid recomputing
    Qp:=Parent(P1`x);
    Qpa<a>:=PolynomialRing(Qp);
    K<a>:=TotallyRamifiedExtension(Qp,a^delta-p);
    format:=recformat<x,b,inf,xt,bt,index>;
    S1:=rec<format|>;                                                    
    S1`inf:=P1`inf;
    S1`x:=Evaluate(xt,a);
    S1`b:=[Evaluate(bt[i],a):i in [1..d]];
  else
    S1:=P1;
  end if;

  if is_bad(P2,data) then
    xt,bt,index:=local_coord(P2,prec,data);
    P2`xt:=xt;       // to avoid recomputing
    P2`bt:=bt;       // to avoid recomputing
    P2`index:=index; // to avoid recomputing
    if not is_bad(P1,data) then
      Qp:=Parent(P2`x);
      Qpa<a>:=PolynomialRing(Qp);
      K<a>:=TotallyRamifiedExtension(Qp,a^delta-p);
    end if;
    format:=recformat<x,b,inf,xt,bt,index>;
    S2:=rec<format|>;                                                    
    S2`inf:=P2`inf;
    S2`x:=Evaluate(xt,a);
    S2`b:=[Evaluate(bt[i],a):i in [1..d]];
  else
    S2:=P2;
  end if;

  tinyP1toS1,NP1toS1:=tiny_integrals_on_basis(P1,S1,data);
  tinyP2toS2,NP2toS2:=tiny_integrals_on_basis(P2,S2,data);

  FS1:=frobenius_pt(S1,data);
  FS2:=frobenius_pt(S2,data);

  tinyS1toFS1,NS1toFS1:=tiny_integrals_on_basis(S1,FS1,data:P:=P1); 
  tinyFS2toS2,NFS2toS2:=tiny_integrals_on_basis(FS2,S2,data:P:=P2); 

  NIP1P2:=Minimum([NP1toS1,NP2toS2,NS1toFS1,NFS2toS2]);

  I:=[];
  for i:=1 to #basis do
    f0iS1,Nf0iS1:=evalf0(f0list[i],S1,data);
    f0iS2,Nf0iS2:=evalf0(f0list[i],S2,data);
    finfiS1,NfinfiS1:=evalfinf(finflist[i],S1,data);
    finfiS2,NfinfiS2:=evalfinf(finflist[i],S2,data);
    fendiS1,NfendiS1:=evalfend(fendlist[i],S1,data);
    fendiS2,NfendiS2:=evalfend(fendlist[i],S2,data);
    NIP1P2:=Minimum([NIP1P2,Nf0iS1,Nf0iS2,NfinfiS1,NfinfiS2,NfendiS1,NfendiS2]);
    I[i]:=(K!f0iS1)-(K!f0iS2)+(K!finfiS1)-(K!finfiS2)+(K!fendiS1)-(K!fendiS2)-(K!tinyS1toFS1[i]+K!tinyFS2toS2[i]);
  end for; 

  mat:=(F-IdentityMatrix(RationalField(),#basis));
  
  mat:=mat^-1;
  valmat:=Minimum([Valuation(e,p):e in Eltseq(mat)]);
  NIP1P2:=NIP1P2+valmat;                            // account for loss of precision multiplying by mat, TODO take error in mat into account as well
  IS1S2:=Vector(I)*Transpose(ChangeRing(mat,K)); 
  
  IP1P2:=IS1S2+ChangeRing(tinyP1toS1,K)-ChangeRing(tinyP2toS2,K);
  IP1P2,Nround:=round_to_Qp(IP1P2);

  assert Nround ge NIP1P2;                          // check that rounding error is within error bound

  NIP1P2:=Ceiling(NIP1P2);

  return IP1P2,NIP1P2;//,tinyP1toS1,NP1toS1,tinyP2toS2,NP2toS2,tinyS1toFS1,NS1toFS1,tinyFS2toS2,NFS2toS2,f0iS1,Nf0iS1,f0iS2,Nf0iS2,finfiS1,NfinfiS1,finfiS2,NfinfiS2,fendiS1,NfendiS1,fendiS2,NfendiS2;
end function;
