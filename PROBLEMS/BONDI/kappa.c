//ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)
//{
  ldouble  kff=kappaCGS2GU(1.7e-25/massGU2CGS(M_PROTON)/massGU2CGS(M_PROTON)* powl(T,-3.5)*rhoGU2CGS(rho));
//relativistic correction (Rybicki & Lightman (5.25))
kff*=(1.+4.4e-10*T);

kff*=rho;  
return kff;

