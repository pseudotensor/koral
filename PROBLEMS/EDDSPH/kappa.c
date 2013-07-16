//ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)

#ifdef OPTTHIN
return 0.;
#endif

#ifdef OPTTHICK
ldouble rhocgs=rhoGU2CGS(rho);
ldouble Tcgs=tempGU2CGS(T);

ldouble ZZsun=1.;
ldouble kappaffcgs=6.4e22*rhocgs/Tcgs/Tcgs/Tcgs/sqrt(Tcgs);
ldouble kappabfcgs=4.8e-24/1.67262158e-24/1.67262158e-24*rhocgs/Tcgs/Tcgs/Tcgs/sqrt(Tcgs)*ZZsun;

return kappaCGS2GU(kappaffcgs)*rho+kappaCGS2GU(kappabfcgs)*rho;
#endif

