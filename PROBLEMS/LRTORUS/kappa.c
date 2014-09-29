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
*kappagasAbs=kappaCGS2GU((6.6e-24/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas)*log(1.+1.6*zeta))*rho;
*kapparadAbs=*kappagasAbs/zeta/zeta/zeta;
	
//Rosseland mean							  
*kappagasRos=kappaCGS2GU((2.1e-25/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas)*(1.-exp(-6.94*zeta)))*rho;
*kapparadRos=*kappagasAbs/zeta/zeta/zeta;

//default
kappa=*kappagasRos;


/*
#ifdef OPACBELLLIN
ldouble X=0.75;
ldouble opaces=0.2*(1.0+X);
ldouble opactot=opacity_BellLin(rhocgs,Tgas);
ldouble opacabs=opactot-opaces;
kappa= kappaCGS2GU(opacabs)*rho;
#else
*/
//ldouble kappaffcgs=6.4e22*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas);
