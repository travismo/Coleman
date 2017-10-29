Fp_points:=function(data);

  // Finds all points on the reduction mod p of the curve given by data

  Q:=data`Q; p:=data`p; d:=Degree(Q); W0:=data`W0; Winf:=data`Winf; 

  Fp:=FiniteField(p); Fpx:=RationalFunctionField(Fp); Fpxy:=PolynomialRing(Fpx);
  f:=Fpxy!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      f:=f+(Fp!Coefficient(Coefficient(Q,i),j))*Fpxy.1^i*Fpx.1^j;
    end for;
  end for;  
  FFp:=FunctionField(f); // function field of curve mod p

  places:=Places(FFp,1);

  b0modp:=[]; // elements of b^0 mod p
  for i:=1 to d do
    f:=FFp!0;
    for j:=1 to d do
      f:=f+(Fpx!W0[i,j])*FFp.1^(j-1);
    end for;
    b0modp[i]:=f;
  end for;

  binfmodp:=[]; // elements of b^inf mod p
  for i:=1 to d do
    f:=FFp!0;
    for j:=1 to d do
      f:=f+(Fpx!Winf[i,j])*FFp.1^(j-1);
    end for;
    binfmodp[i]:=f;
  end for;

  Fppts:=[];
  for i:=1 to #places do
    if Valuation(FFp!Fpx.1,places[i]) ge 0 then
      if Valuation(FFp!Fpx.1-Evaluate(FFp!Fpx.1,places[i]),places[i]) eq 1 then
        index:=0;
      else
        j:=1;
        done:=false;
        while not done and j le d do
          if Valuation(b0modp[j]-Evaluate(b0modp[j],places[i]),places[i]) eq 1 then
            done:=true;
            index:=j;
          end if;
          j:=j+1;
        end while;
      end if;
      Fppts:=Append(Fppts,[*false,Evaluate(FFp!Fpx.1,places[i]),[Fp!Evaluate(b0modp[j],places[i]): j in [1..d]],index*]);    
    else
      if Valuation(FFp!(1/Fpx.1)-Evaluate(FFp!(1/Fpx.1),places[i]),places[i]) eq 1 then
        index:=0;
      else
        j:=1;
        done:=false;
        while not done and j le d do
          if Valuation(binfmodp[j]-Evaluate(binfmodp[j],places[i]),places[i]) eq 1 then
            done:=true;
            index:=j;
          end if;
          j:=j+1;
        end while;
      end if;
      Fppts:=Append(Fppts,[*true,Evaluate(FFp!(1/Fpx.1),places[i]),[Fp!Evaluate(binfmodp[j],places[i]): j in [1..d]],index*]);
    end if;
  end for;

  return Fppts;

end function;


Qp_points:=function(data:points:=[]);

  // For every point on the reduction mod p of the curve given by data,
  // a Qp point on the curve is returned that reduces to it. Optionally,
  // an (incomplete) list pts can be specified by the user which will
  // then be completed.
  
  Q:=data`Q; p:=data`p; N:=data`N; W0:=data`W0; Winf:=data`Winf; 
  d:=Degree(Q); Fp:=FiniteField(p); Qp:=pAdicField(p,N); Qpy:=PolynomialRing(Qp);
  Qx:=RationalFunctionField(RationalField()); Qxy:=PolynomialRing(Qx);

  Fppts:=Fp_points(data);
  Qppts:=[];

  f:=Qxy!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      f:=f+Coefficient(Coefficient(Q,i),j)*Qxy.1^i*Qx.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  b0fun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+W0[i,j]*FF.1^(j-1);
    end for;
    b0fun[i]:=bi;
  end for;

  binffun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+Winf[i,j]*FF.1^(j-1);
    end for;
    binffun[i]:=bi;
  end for;

  for i:=1 to #Fppts do
    
    Fppoint:=Fppts[i];
    j:=1;
    done:=false;
    while not done and j le #points do
      if (Fppoint[1] eq points[j]`inf) and (Fp!(points[j]`x)-Fppoint[2] eq 0) and ([Fp!(points[j]`b)[k]:k in [1..d]] eq Fppoint[3]) then
        done:=true;
        P:=points[j];
      end if;
      j:=j+1; 
    end while;
    
    if not done then
      
      if Fppoint[1] then // infinite point
        
        inf:=true;
        
        if Fppoint[4] eq 0 then // x - point[1] local coordinate
          x:=Qp!Fppoint[2];
          b:=[];
          for j:=1 to d do
            bj:=binffun[j];
            poly:=minpoly(FF!(1/Qx.1),bj);
            C:=Coefficients(poly);
            D:=[];
            for k:=1 to #C do
              D[k]:=Evaluate(C[k],x); 
            end for;
            fy:=Qpy!D;
            zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
            done:=false;
            k:=1;
            while not done and k le #zeros do
              if (Fp!zeros[k][1]-Fppoint[3][j] eq 0) then 
                done:=true;
                b[j]:=zeros[k][1];
              end if;
              k:=k+1;
            end while;
          end for;
        else // x-point[1] not local coordinate
          index:=Fppoint[4];
          bindex:=Qp!Fppoint[3][index];
          poly:=minpoly(binffun[index],FF!(1/Qx.1));
          C:=Coefficients(poly);
          D:=[];
          for k:=1 to #C do
            D[k]:=Evaluate(C[k],bindex); 
          end for;
          fy:=Qpy!D;
          zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
          done:=false;
          k:=1;
          while not done and k le #zeros do
            if (Fp!zeros[k][1]-Fppoint[2] eq 0) then 
              done:=true;
              x:=zeros[k][1];
            end if;
            k:=k+1;
          end while;
          b:=[];
          for j:=1 to d do
            if j eq index then
              b[j]:=bindex;
            else
              poly:=minpoly(binffun[index],binffun[j]);
              C:=Coefficients(poly);
              D:=[];
              for k:=1 to #C do
                D[k]:=Evaluate(C[k],bindex); 
              end for;
              fy:=Qpy!D;
              zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
              done:=false;
              k:=1;
              while not done and k le #zeros do
                if (Fp!zeros[k][1]-Fppoint[3][j] eq 0) then 
                  done:=true;
                  b[j]:=zeros[k][1];
                end if;
                k:=k+1;
              end while;
            end if;
          end for; 
        end if;
      else // finite point
        inf:=false;
        if Fppoint[4] eq 0 then // x - point[1] local coordinate
          x:=Qp!Fppoint[2];
          b:=[];
          for j:=1 to d do
            bj:=b0fun[j];
            poly:=minpoly(FF!Qx.1,bj);
            C:=Coefficients(poly);
            D:=[];
            for k:=1 to #C do
              D[k]:=Evaluate(C[k],x); 
            end for;
            fy:=Qpy!D;
            zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
            done:=false;
            k:=1;
            while not done and k le #zeros do
              if (Fp!zeros[k][1]-Fppoint[3][j] eq 0) then 
                done:=true;
                b[j]:=zeros[k][1];
              end if;
              k:=k+1;
            end while;
          end for;
        else // x-point[1] not local coordinate
          index:=Fppoint[4];
          bindex:=Qp!Fppoint[3][index];
          poly:=minpoly(b0fun[index],FF!Qx.1);
          C:=Coefficients(poly);
          D:=[];
          for k:=1 to #C do
            D[k]:=Evaluate(C[k],bindex); 
          end for;
          fy:=Qpy!D;
          zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
          done:=false;
          k:=1;
          while not done and k le #zeros do
            if (Fp!zeros[k][1]-Fppoint[2] eq 0) then 
              done:=true;
              x:=zeros[k][1];
            end if;
            k:=k+1;
          end while;
          b:=[];
          for j:=1 to d do
            if j eq index then
              b[j]:=bindex;
            else
              poly:=minpoly(b0fun[index],b0fun[j]);
              C:=Coefficients(poly);
              D:=[];
              for k:=1 to #C do
                D[k]:=Evaluate(C[k],bindex); 
              end for;
              fy:=Qpy!D;
              zeros:=Roots(fy); // Hensel lifting gives problems here, since Hensel condition not always satisfied
              done:=false;
              k:=1;
              while not done and k le #zeros do
                if (Fp!zeros[k][1]-Fppoint[3][j] eq 0) then 
                  done:=true;
                  b[j]:=zeros[k][1];
                end if;
                k:=k+1;
              end while;
            end if;
          end for;
        end if;
        
      end if;
    
    P:=set_bad_point(x,b,inf,data);

    end if;
    
    Qppts:=Append(Qppts,P);
  
  end for;

  return Qppts;

end function;


Q_points:=function(data,bound);

  // Returns a list (not guaranteed to be complete) of Q-rational points 
  // upto height bound on the curve given by data.

  Q:=data`Q; p:=data`p; N:=data`N; r:=data`r; W0:=data`W0; Winf:=data`Winf;
  d:=Degree(Q);

  pointlist:=[];

  A2:=AffineSpace(RationalField(),2);
  Qxy:=PolynomialRing(RationalField(),2);
  
  QA2:=Qxy!0;
  C:=Coefficients(Q);
  for i:=1 to #C do
    D:=Coefficients(C[i]);
    for j:=1 to #D do
      QA2:=QA2+D[j]*(Qxy.1)^(j-1)*(Qxy.2)^(i-1);
    end for;
  end for;
  
  X:=Scheme(A2,QA2);
  pts:=PointSearch(X,bound);
  xvalues:=[];
  for i:=1 to #pts do
    if not pts[i][1] in xvalues then
      xvalues:=Append(xvalues,pts[i][1]);
    end if;
  end for;

  Qx:=RationalFunctionField(RationalField());
  Qxy:=PolynomialRing(Qx);
  Qfun:=Qxy!0;
  C:=Coefficients(Q);
  for i:=1 to #C do
    D:=Coefficients(C[i]);
    for j:=1 to #D do
      Qfun:=Qfun+D[j]*(Qx.1)^(j-1)*(Qxy.1)^(i-1);
    end for;
  end for;
  FF:=FunctionField(Qfun);

  b0fun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+W0[i,j]*FF.1^(j-1);
    end for;
    b0fun[i]:=bi;
  end for;

  binffun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+Winf[i,j]*FF.1^(j-1);
    end for;
    binffun[i]:=bi;
  end for;

  for i:=1 to #xvalues do
    places:=Decomposition(FF,Zeros(Qx.1-xvalues[i])[1]);
    if Valuation(xvalues[i],p) ge 0 then
      for j:=1 to #places do
        if Degree(places[j]) eq 1 then
          x:=xvalues[i];
          b:=[];
          for k:=1 to d do
            b[k]:=Evaluate(b0fun[k],places[j]);
          end for;
          P:=set_bad_point(x,b,false,data);
          pointlist:=Append(pointlist,P);
        end if;
      end for;
    else
      for j:=1 to #places do
        if Degree(places[j]) eq 1 then
          x:=1/xvalues[i];
          b:=[];
          for k:=1 to d do
            b[k]:=Evaluate(binffun[k],places[j]);
          end for;
          P:=set_bad_point(x,b,true,data);
          pointlist:=Append(pointlist,P);
        end if;
      end for;
    end if;
  end for; 

  places:=InfinitePlaces(FF);
  for i:=1 to #places do
    if Degree(places[i]) eq 1 then
      x:=0;
      b:=[];
      for j:=1 to d do
        b[j]:=Evaluate(binffun[j],places[i]);
      end for;
      P:=set_bad_point(x,b,true,data);
      pointlist:=Append(pointlist,P);
    end if;
  end for; 

  return pointlist;

end function;


my_roots_Zpt:=function(f)

  // Custom function to compute the roots of a polynomial
  // f over Z_p since the Magma intrinsic requires the leading
  // coefficient to be a unit (which is usually not the case 
  // for us).

  if f eq 0 then
    error "Polynomial has to be non-zero";
  end if;

  Zps:=Parent(f);
  Zp:=BaseRing(Zps);
  Fp:=ResidueClassField(Zp);
  Fps:=PolynomialRing(Fp);
  p:=Characteristic(Fp);
  
  Nf:=Precision(Zp);
  val:=Minimum([Valuation(e):e in Eltseq(f)]);

  Zp:=ChangePrecision(Zp,Nf-val);
  Zps:=PolynomialRing(Zp);
  
  f:=Zps![e/p^val :e in Eltseq(f)];
  i:=0;
  done:=false;
  while not done do
    if Coefficient(f,i) ne 0 then
      done:=true;
    end if;
    i:=i+1;
  end while;
  if i gt 0 then
    coefs:=Coefficients(f);
    for j:=i-1 to 1 by -1 do
      coefs:=Remove(coefs,1);
    end for;
    f:=Zps!coefs;
    zero:=true;
  end if;

  modproots:=Roots(Fps!f);
  Fproots:=[];
  for i:=1 to #modproots do
    Fproots[i]:=modproots[i][1];
  end for;
  Zproots:=[[*Zp!e,1*]:e in Fproots];

  i:=1;
  while i le #Zproots do
    z:=Zproots[i][1];
    Nz:=Zproots[i][2];
    v1:=Valuation(Evaluate(f,z));
    v2:=Valuation(Evaluate(Derivative(f),z)); 
    if (not v1 gt 2*v2) and (v1 lt Nf-val) then
      Zproots:=Remove(Zproots,i);
      znew:=z+p^Nz*Zps.1;
      g:=Fps![e/p^(Nz): e in Coefficients(Evaluate(f,znew))];
      if g ne 0 then
        Fproots:=Roots(g);
      else
        Fproots:=[[e,1]: e in Fp];
      end if;
      for j:=1 to #Fproots do
        Zproots:=Insert(Zproots,i,[*z+p^Nz*(Zp!Fproots[j][1]),Nz+1*]);
      end for;
    else
      i:=i+1;
    end if;
  end while;

  for i:=1 to #Zproots do
    z:=Zproots[i][1];
    Nz:=Zproots[i][2];
    v1:=Valuation(Evaluate(f,z));
    v2:=Valuation(Evaluate(Derivative(f),z));
    if (v1 lt Nf-val) then
      z:=HenselLift(f,z);
      Zproots[i][1]:=z;
      Zproots[i][2]:=Nf-val-v2;
    else
      Zproots[i][2]:=Nf-val-v2;
    end if;
  end for;

  if zero then
    Zproots:=Append(Zproots,[*Zp!0,Nf-val*]);
  end if;

  return Zproots;

end function;


basis_kernel:=function(A)

  // Compute a basis for the kernel of the matrix A over Qp

  val:=Minimum([0] cat [Valuation(x) : x in Eltseq(A)]); 
  Qp:=BaseRing(A);
  N:=Precision(Qp);
  p:=Prime(Qp);

  A:=p^(-val)*A; 
  N:=N-val;

  Zp:=pAdicRing(p,N);
  row:=NumberOfRows(A);
  col:=NumberOfColumns(A);
  matpN:=ZeroMatrix(Zp,row,col);
  for i:=1 to row do
    for j:=1 to col do
      matpN[i,j]:=Zp!A[i,j];
    end for;
  end for;
  S,P1,P2:=SmithForm(matpN);            
 
  b:=[];
  for i:=Rank(S)+1 to row do
    Append(~b,P1[i]);
  end for;
  if #b gt 0 then
    b:=RowSequence(HermiteForm(Matrix(b)));
  end if;
 
  b:=RowSequence(ChangeRing(Matrix(b),Qp));

  return b;
end function;


vanishing_differentials:=function(points,data:e:=1);

  // Compute the regular one forms of which the 
  // integrals vanish between all points in points.

  Q:=data`Q; p:=data`p;
  
  g:=genus(Q,p);

  IP1Pi:=[];
  NIP1Pi:=[];
  for i:=1 to #points do
    Ci,Ni:=coleman_integrals_on_basis(points[1],points[i],data:e:=e);
    IP1Pi[i]:=Ci;
    NIP1Pi[i]:=Ni;
  end for;

  Nint:=Minimum(NIP1Pi);

  K:=pAdicField(p,Nint);
  M:=ZeroMatrix(K,g,#points);
  for i:=1 to g do
    for j:=1 to #points do
      M[i,j]:=K!reduce_mod_pN_Q(Rationals()!IP1Pi[j][i],p,Nint);
    end for;
  end for;

  v:=basis_kernel(M);

  return v,Nint;

end function;


zeros_on_disk:=function(P1,P2,v,data:prec:=0,e:=1);

  // Find all common zeros of the integrals of the v[i] (vectors 
  // of length g) from P1 to points in the residue disk of P2.

  Q:=data`Q; p:=data`p; N:=data`N; 

  g:=genus(Q,p);

  IP1P2,NIP1P2:=coleman_integrals_on_basis(P1,P2,data:e:=e);
  tinyP2toz,xt,bt,NP2toz:=tiny_integrals_on_basis_to_z(P2,data:prec:=prec);

  Zp:=pAdicRing(p,NIP1P2);
  Zpt:=PolynomialRing(Zp);

  zerolist:=[];
  for i:=1 to #v do
    f:=Parent(tinyP2toz[1])!0;
    for j:=1 to g do
      f:=f+v[i][j]*(IP1P2[j]+tinyP2toz[j]);
    end for;
    h:=Zpt!0;
    for j:=0 to Degree(f) do
      h:=h+IntegerRing()!(p^j*(RationalField()!Coefficient(f,j)))*Zpt.1^j;
    end for; 
    zeros:=my_roots_Zpt(h);
    zerolist:=Append(zerolist,zeros);
  end for;

  zeroseq:=[];
  for i:=1 to #zerolist[1] do
    allzero:=true;
    for j:=2 to #zerolist do
      found:=false;
      for k:=1 to #zerolist[j] do
        if Valuation(zerolist[j][k][1]-zerolist[1][i][1]) ge Minimum(zerolist[j][k][2],zerolist[1][i][2]) then
          found:=true;
        end if;
      end for;
      if not found then
        allzero:=false;
      end if;
    end for;
    if allzero then
      zeroseq:=Append(zeroseq,zerolist[1][i][1]);
    end if; 
  end for;

  pointlist:=[];
  for i:=1 to #zeroseq do
    z:=zeroseq[i];
    x:=Evaluate(xt,p*z); 
    b:=Eltseq(Evaluate(bt,p*z));
    inf:=P2`inf; 
    P:=set_bad_point(x,b,P2`inf,data);
    pointlist:=Append(pointlist,P);
  end for;

  return pointlist;

end function;
