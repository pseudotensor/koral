return 0.;

return 10000.*.1*rho;

/*
ldouble rhocgs=rhoGU2CGS(rho);
ldouble Tcgs=tempGU2CGS(T);

ldouble ZZsun=1.;
ldouble kappaffcgs=6.4e22*rho/Tcgs/Tcgs/Tcgs/sqrt(Tcgs);
ldouble kappabfcgs=4.8e-24/1.67262158e-24/1.67262158e-24*rho/Tcgs/Tcgs/Tcgs/sqrt(Tcgs)*ZZsun;

return kappaCGS2GU(kappaffcgs)*rho+kappaCGS2GU(kappabfcgs)*rho;
*/
