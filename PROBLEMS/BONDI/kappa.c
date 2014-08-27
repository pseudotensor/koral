//absorption opacities

ldouble rhocgs=rhoGU2CGS(rho);
ldouble Tgas=T;
ldouble Trad=Tgas;

#ifdef NCOMPTONIZATION //adjust the Rosseland mean for the color temperature of radiation
Trad = calc_ncompt_Thatrad_full(pp,ggg);
#endif

ldouble zeta = Trad/Tgas;


		   
//absorbtion mean
ldouble mpcgs=1.67262158e-24;
*kappagasAbs=kappaCGS2GU((6.6e-24/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas)*log(1.+1.6*zeta))*rho*(1.+4.4e-10*Tgas);
*kapparadAbs=*kappagasAbs/zeta/zeta/zeta;
	
//Roseland mean							  
*kappagasRos=kappaCGS2GU((2.1e-25/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas)*(1.-exp(-6.94*zeta)))*rho*(1.+4.4e-10*Tgas);
*kapparadRos=*kappagasAbs/zeta/zeta/zeta;

//default
kappa=*kappagasRos;

