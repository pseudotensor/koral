//absorption

ldouble rhocgs=rhoGU2CGS(rho);
ldouble Tcgs=tempGU2CGS(T);

#ifdef OPACBELLLIN

ldouble X=0.75;
ldouble opaces=0.2*(1.0+X);
ldouble opactot=opacity_BellLin(rhocgs,Tcgs);
ldouble opacabs=opactot-opaces;
return kappaCGS2GU(opacabs)*rho;

#else

ldouble kappaffcgs=6.4e22*rhocgs/Tcgs/Tcgs/Tcgs/sqrt(Tcgs);
return kappaCGS2GU(kappaffcgs)*rho;

#endif


