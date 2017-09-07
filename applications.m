my_roots_Zpt:=function(f)

  // Custom function to compute the roots of a polynomial
  // f over Z_p since the Magma intrinsic requires the leading
  // coefficient to be a unit (which is usually not the case 
  // for us).

  val:=Minimum([Valuation(e):e in Eltseq(f)]);

  Zps<s>:=Parent(f);
  Zp:=BaseRing(Zps);
  Fp:=ResidueClassField(Zp);
  Fps:=PolynomialRing(Fp);
  p:=Characteristic(Fp);
  
  Nf:=Precision(Zp);
  Zp:=ChangePrecision(Zp,Nf-val);
  Zps:=PolynomialRing(Zp);
  
  f:=Zps![e/p^val :e in Eltseq(f)];
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
    end if;
  end for;

  return Zproots;

end function;


torsion_points_in_disk:=function(P1,P2,data:prec:=0,delta:=1)

  // Finds the points P in residue disk of P2 such that the integral of
  // all regular 1-forms from P1 to P vanishes. Note that this set of points
  // contains the set of points P for which P-P1 is torsion in the Jacobian.

  p:=data`p; N:=data`N; basis:=data`basis; g:=IntegerRing()!(#Eltseq(basis)/2);

  IP1P2,NIP1P2:=coleman_integrals_on_basis(P1,P2,data:delta:=delta);
  tinyP2toz,xt,bt,NP2toz:=tiny_integrals_on_basis_to_z(P2,data:prec:=prec);

  IP1toz:=IP1P2+tinyP2toz;
  NIP1toz:=Minimum(NIP1P2,NP2toz);
  
  prec:=Maximum([Degree(e): e in Eltseq(tinyP2toz)]);
  Zp:=pAdicRing(p,NIP1P2); 
  Zpt:=PowerSeriesRing(Zp,prec);
  Zps:=PolynomialRing(Zp);

  zerolist:=[];
  for i:=1 to g do
    f:=Eltseq(IP1toz)[i];
    g:=Zpt!0;
    for j:=0 to Degree(IP1toz[i]) do
      if Valuation(Coefficient(f,j)) lt NIP1P2 then
        g:=g+Zp!IntegerRing()!(p^j*(RationalField()!Coefficient(f,j)))*Zpt.1^j;
      end if;    
    end for;
    h:=(Zps!Coefficients(g));
    zeros:=my_roots_Zpt(h);
    val:=Valuation(g);
    if val gt 0 then
      zeros:=Append(zeros,[*Zp!0,Ceiling(NIP1P2/val)*]);
    end if;
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

  format:=recformat<x,b,inf,xt,bt,index>;
  pointlist:=[];
  for i:=1 to #zeroseq do
    z:=zeroseq[i];
    P:=rec<format|>;
    P`inf:=P2`inf;
    P`x:=Evaluate(xt,p*z);
    P`b:=Vector(Evaluate(bt,p*z));
    pointlist:=Append(pointlist,P);
  end for;
  
  return pointlist,NIP1toz;

end function;
