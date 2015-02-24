//absorption opacity
//free-free + bound-free

ldouble rhocgs=rhoGU2CGS(rho);
ldouble Tcgs=tempGU2CGS(T);

ldouble ZZsun=1.;
ldouble kappaffcgs=6.4e22*rho/Tcgs/Tcgs/Tcgs/sqrt(Tcgs);
//ldouble kappabfcgs=4.8e-24/1.67262158e-24/1.67262158e-24*rho/Tcgs/Tcgs/Tcgs/sqrt(Tcgs)*ZZsun;

kappa=1.e-1;

//does not work with FIXALLBUTEMP
//kappa*=step_function(log(rho)-(-14.),1.);


return kappaCGS2GU(kappa)*rho;
//return kappaCGS2GU(kappaffcgs)*rho;

