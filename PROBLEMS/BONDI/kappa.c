//absorption opacities

ldouble rhocgs=rhoGU2CGS(rho);
ldouble Tgas=T;
ldouble Trad=Tgas;

#ifdef NCOMPTONIZATION //the color temperature of radiation
Trad = calc_ncompt_Thatrad_full(pp,ggg);
#endif

ldouble zeta = Trad/Tgas;
//test
zeta=1.;

ldouble mpcgs=1.67262158e-24;

#ifndef SKIPFANCYOPACITIES		   

//absorbtion mean
ldouble kappaff,kappabe;
kappaff=kappaCGS2GU((6.6e-24/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas)*log(1.+1.6*zeta))*rho*(1.+4.4e-10*Tgas);
kappabe=kappaCGS2GU((8.3e-15/mpcgs/mpcgs)*rhocgs*pow(Tgas,-1.7)/Tgas/Tgas/Tgas*log(1.+1.6*zeta))*rho*(1.+4.4e-10*Tgas);

*kappagasAbs=kappaff+kappabe;
*kapparadAbs=*kappagasAbs/zeta/zeta/zeta;
	
//Roseland mean - not used at all							  
*kappagasRos=kappaCGS2GU((2.1e-25/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas)*(1.-exp(-6.94*zeta)))*rho*(1.+4.4e-10*Tgas);
*kapparadRos=*kappagasAbs/zeta/zeta/zeta;

//default
kappa=*kappagasAbs;

//test
*kappagasAbs=kappaff;
*kapparadAbs=kappabe;

#else
//the simplest
kappa=kappaCGS2GU((6.6e-24/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas))*rho*(1.+4.4e-10*Tgas);



#endif



